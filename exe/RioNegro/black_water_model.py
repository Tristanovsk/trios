import os
import numpy as np
import pandas as pd

import matplotlib as mpl
import matplotlib.pyplot as plt

plt.ioff()
mpl.rcParams.update({'font.size': 18})
import statsmodels.api as sm
from sklearn.linear_model import LinearRegression
from statsmodels.sandbox.regression.predstd import wls_prediction_std
from scipy import odr as odr
from scipy import stats as scistats

import RTxploitation as rt
iopw = rt.auxdata.iopw().get_iopw

opj= os.path.join
idir = os.path.abspath('/DATA/OBS2CO/data/rogerio')
figdir = opj(idir,'fig')
file = opj(idir,'data/Simulation_Rrs_OSOAA_TSS_COD_aCDOM.xlsx')

data = pd.read_excel(file)
data.sort_values(by="TSS (mg L)", inplace=True)
data = data.set_index(['Station', 'Date'])

acdom = data.iloc[:, 17]
acdom_sd = data.iloc[:, 18]
#remove sparse data of CDOM
data = data.iloc[:,:-2]
params=['SPM','DOC']
param=params[0]
if param == 'DOC':
    data.dropna(inplace=True)
month = data.index.get_level_values(1).month
spm = data.iloc[:, 0]  # .values.reshape(-1, 1)
spm_sd = data.iloc[:, 1]
doc = data.iloc[:, 15]
doc_sd = data.iloc[:, 16]

def Rrs2rrs(Rrs):
        return Rrs / (0.52 + 1.7 * Rrs)

def inv_gordon88(rrs,g0=0.089,g1=0.125):
        deltas = g0**2 + 4 * g1 * rrs
        return (-g0 + np.sqrt(deltas)) / (2*g1)

def black_water_model_old(B, x , wl=550):
    aw,bbw = iopw(wl)
    rrs = Rrs2rrs(x)
    u = inv_gordon88(rrs)
    return (B[0])**2 * (u - np.abs(B[1])  )

def black_water_model(B, x, wl=550):

    rrs = Rrs2rrs(x)
    u = inv_gordon88(rrs)
    aw,bbw = iopw(wl)
    N = (u*(aw+bbw+B[0])-bbw)/ (B[1]- u *(B[1]+B[2]))

    return B[3] * N

def linear(B, x):
    '''Linear function y = m*x + b'''
    return B[0] + B[1] * x

def exponential(B, x):
    #return B[0] - (B[0]-B[1]) * np.exp(-B[2] * x)
    return B[0] * np.exp(B[1] * x) + B[2]

def hyperbolic(B,x,b=1):
    return B[1]/(1+B[0]*x)**(1/b)

def poly(B, x):
    return B[0] + B[1] * x + B[2] * x ** 2  # + B[3]*x**3

def confidence_interval(xn, model, res, nstd=1):
    '''

    :param xn: x-axis data for computation
    :param model: numerical model used for the fit
    :param res: output from scipy.odr.run
    :param nstd: number of sigma to compute confidence interval
    :return: data up and data down
    '''
    '''
    
    :param res: output from scipy.odr.run
    :param nstd: number of sigma to compute confidence interval
    :return: 
    '''

    popt_up = res.beta + nstd * res.sd_beta
    popt_dw = res.beta - nstd * res.sd_beta

    return model(popt_up, xn), model(popt_dw, xn)

if param == 'SPM':
    y=spm
    y_sd=spm_sd
    fit_eq='$g_0={:.1f};\ g_1={:.1f};\ C={:.1f}$'
    ylabel='SPM (mg/L)'
    # gordon values
    g0 = 0.089
    g1 = 0.125
    model,beta0=black_water_model_old,[150,10]
    #model,beta0=black_water_model,[1,1,1,10]
else:
    y=doc
    y_sd=doc_sd
    fit_eq='y = {:.1f}/(1+{:.1f}$x)$\n res_var = {:.1f}'
    ylabel='DOC (mg/L)'
    model,beta0=hyperbolic,[10,10]

c=doc
cmap = plt.cm.get_cmap("Spectral").reversed()
stats = pd.DataFrame(columns=['band','algo','b0','b1','b2','sig_b0','sig_b1','sig_b2'])


# plot Rrs = f(SPM)
wls = [560, 665, 704,740, 782, 864]
fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(20, 12))
for i, ax in enumerate(axs.ravel()):
    wl=wls[i]

    Rrs = data.iloc[:, 2 * i + 2]
    print(Rrs.name,wl)
    Rrs_sd = data.iloc[:, 2 * i + 3]

    x,x_sd = Rrs, Rrs_sd

    ax.errorbar(x, y, xerr=x_sd, yerr=y_sd, fmt='o', color='grey', ecolor='black', alpha=0.8)
    # to put color in scatter dots
    #ax.scatter(x, y, c=c, s=55, cmap=cmap)
    ax.set_title(x.name)
    ax.set_xlabel(r'$R_{rs}\ (sr^{-1})$')
    ax.set_ylabel(ylabel)

    testdata = odr.RealData(x, y, sx=x_sd, sy=y_sd)
    _odr = odr.ODR(testdata, odr.Model(model, extra_args=[wl]), beta0=beta0) #

    # fit with ordinary least square (OLS, fit_type=2 )
    _odr.set_job(fit_type=2)
    res = _odr.run()

    res.pprint()
    xn = np.linspace(0, np.max(x) * 1.25, 500)
    yn = model(res.beta, xn)
    fit_up, fit_dw = confidence_interval(xn, model, res)

    #stats.loc[2*i]=np.concatenate([[Rrs.name, 'OLS'], res.beta, res.sd_beta])


    ax.plot(xn, yn, 'b--', label='OLS',#; '+fit_eq.format(*res.beta),
            linewidth=2)
    ax.fill_between(xn, fit_up, fit_dw, alpha=.25, facecolor="b")  # ,label="1-sigma interval")

    # fit with Orthogonal distance regression (ODR, fit_type = 0)
    _odr.set_job(fit_type=0)
    res = _odr.run()
    res.pprint()
    xn = np.linspace(0, np.max(x) * 1.25, 500)
    yn = model(res.beta, xn)
    fit_up, fit_dw = confidence_interval(xn, model, res)
    print()
    #stats.loc[2*i+1]=np.concatenate([[Rrs.name, 'ODR'], res.beta, res.sd_beta])

    ax.plot(xn, yn, 'r-', label='ODR',#; '+fit_eq.format(*res.beta),
            linewidth=2)
    ax.fill_between(xn, fit_up, fit_dw, alpha=.25, facecolor="r")  # ,label="1-sigma interval")
    ax.legend()
    ax.set_ylim([0,40])
plt.suptitle('Fit based on modified Gordon model')
stats.to_csv(opj(idir,'stats_SPM_Rrs.csv'), float_format='%.2f',index=False)
plt.tight_layout(rect=[0.0, 0.0, 0.99, 0.94])
fig.savefig(opj(figdir,'blackwater_'+param+'_vs_Rrs.png'), dpi=300)
plt.show()

plt.close()
stats= stats.set_index(['band','algo'])
stats = stats.astype('float')

#------------------------------------------------
# Error histograms
#------------------------------------------------

fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(20, 12))
axs=axs.ravel()
i=0
bandwidth=0.4
for g,d in data.loc[:,['Rrs_sim(B3)','Rrs_sim(B4)','Rrs_sim(B5)','Rrs_sim(B6)','Rrs_sim(B7)','Rrs_sim(B8a)']].iteritems():
    print(g)
    ax=axs[i]
    for method in ('OLS','ODR'):
        c='red'
        if method == 'OLS':
            c='blue'
        stat=stats.loc[(g,method)]

        err = model([stat.b0, stat.b1,stat.b2], d) - y.values.ravel()

        kde = scistats.gaussian_kde(err,bw_method=bandwidth)# / x.std(ddof=1))
        xx = np.linspace(-12,12, 1000)

        ax.hist(err, bins=np.arange(-15,15,bandwidth), color=c,label=method, alpha=0.5, density=True,rwidth=0.8)
        ax.plot(xx, kde.evaluate(xx),color=c,label=method,lw=2)
    ax.axes.get_yaxis().set_visible(False)
    ax.set_title(g)
    ax.set_xlim(-12,12)
    ax.set_xlabel(r'Retrieved - Measured SSSC (mg/L)')
    ax.set_ylabel(r'N')
    ax.legend()
    print(y,)
    i+=1
plt.suptitle('Normalized histogram of errors')
plt.tight_layout(rect=[0.05, 0.05, 0.99, 0.95])
#fig.savefig(opj(figdir,'blackwater_histogram_SPM-err-retrieval.png'), dpi=200)
plt.show()

#------------------------------------------------
# compare with uncertainty
#------------------------------------------------
def sig_param(Rrs,sig_a,sig_b,cov_ab):
    return ((Rrs*sig_a)**2 + sig_b**2 + 2*Rrs*cov_ab)**0.5

fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(20, 12))
axs=axs.ravel()
i=0

for g,d in data.loc[:,['Rrs_sim(B3)','Rrs_sim(B4)','Rrs_sim(B5)','Rrs_sim(B6)','Rrs_sim(B7)','Rrs_sim(B8a)']].iteritems():
    print(g)
    ax=axs[i]
    for method in ('OLS','ODR'):
        c='red'
        if method == 'OLS':
            c='blue'
        stat=stats.loc[(g,method)]

        est = model([stat.b0, stat.b1,stat.b2], d)
        slope, intercept, r_value, p_value, std_err = scistats.linregress(y, est)
        ax.plot([-10, 100], intercept + slope * np.array([-10, 100]), color=c, lw=1.6)
        #sig_est = sig_param(d,stat.sig_a,stat.sig_b,stat.cov_ab)
        #yerr = sig_est,
        ax.errorbar(y, est, xerr=y_sd,  fmt='o', color=c, ecolor=c, label=method+'; y={:.1f}x+{:.1f}; $R^2$={:.3f}'.format(slope, intercept,r_value**2),alpha=0.8)


    ax.set_title(g)

    ax.set_xlabel(r'Measured SSSC (mg/L)')
    ax.set_ylabel(r'Retrieved SSSC (mg/L)')
    ax.plot([-10, 100], [-10, 100], '--', color="grey", lw=1.6)
    ax.set_xlim(0,30)
    ax.set_ylim(0,30)

    ax.legend()
    print(y,)
    i+=1
plt.suptitle('Comparison retrievals')
plt.tight_layout(rect=[0.05, 0.05, 0.99, 0.95])
fig.savefig(opj(figdir,'blackwater_compar_SPM-retrieval.png'), dpi=200)
plt.show()

