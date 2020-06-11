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

def black_water_model(B, x):

    rrs = Rrs2rrs(x)
    return B[2] * (-B[0]+np.sqrt(B[0]**2-4*B[1]*rrs)) / (2*B[1])

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


models = [(hyperbolic,2),(linear, 2),(exponential,3)]  # (poly,3),,


if param == 'SPM':
    y=spm
    y_sd=spm_sd
    fit_eq='y = {:.1f}x + {:.1f}\n res_var = {:.1f}'
    ylabel='SPM (mg/L)'
    model,N=models[1]
else:
    y=doc
    y_sd=doc_sd
    fit_eq='y = {:.1f}/(1+{:.1f}$x)$\n res_var = {:.1f}'
    ylabel='DOC (mg/L)'
    model,N=models[0]

c=doc
cmap = plt.cm.get_cmap("Spectral").reversed()
stats = pd.DataFrame(columns=['band','algo','a','b','rel_sd_a','rel_sd_b','sig_a','sig_b','cov_ab'])


# plot Rrs = f(SPM)

fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(20, 12))
for i, ax in enumerate(axs.ravel()):
    print(i)
    Rrs = data.iloc[:, 2 * i + 2]
    print(Rrs.name)
    Rrs_sd = data.iloc[:, 2 * i + 3]

    x,x_sd = Rrs, Rrs_sd

    ax.errorbar(x, y, xerr=x_sd, yerr=y_sd, fmt='o', color='grey', ecolor='black', alpha=0.8)
    # to put color in scatter dots
    #ax.scatter(x, y, c=c, s=55, cmap=cmap)
    ax.set_title(x.name)
    ax.set_xlabel(r'$R_{rs}\ (sr^{-1})$')
    ax.set_ylabel(ylabel)

    testdata = odr.RealData(x, y, sx=x_sd, sy=y_sd)
    _odr = odr.ODR(testdata, odr.Model(model), beta0=[10.]*N)

    # fit with ordinary least square (OLS, fit_type=2 )
    _odr.set_job(fit_type=2)
    res = _odr.run()

    res.pprint()
    xn = np.linspace(0, np.max(x) * 1.25, 50)
    yn = model(res.beta, xn)
    fit_up, fit_dw = confidence_interval(xn, model, res)
    a, b, sd_a, sd_b, cov_ab = res.beta[1], res.beta[0], res.sd_beta[1], res.sd_beta[0],res.cov_beta[0,1]
    stats.loc[2*i]=[Rrs.name, 'OLS', a,b,100*sd_a/a,100*sd_b/b, sd_a, sd_b, cov_ab]

    ax.plot(xn, yn, 'b--', label='OLS; '+fit_eq.format(a, b, res.res_var),
            linewidth=2)
    #ax.fill_between(xn, fit_up, fit_dw, alpha=.25, facecolor="b")  # ,label="1-sigma interval")

    # fit with Orthogonal distance regression (ODR, fit_type = 0)
    _odr.set_job(fit_type=0)
    res = _odr.run()
    res.pprint()
    xn = np.linspace(0, np.max(x) * 1.25, 50)
    yn = model(res.beta, xn)
    fit_up, fit_dw = confidence_interval(xn, model, res)

    a, b, sd_a, sd_b, cov_ab = res.beta[1], res.beta[0], res.sd_beta[1], res.sd_beta[0],res.cov_beta[0,1]
    stats.loc[2*i+1]=[Rrs.name, 'ODR', a,b,100*sd_a/a,100*sd_b/b, sd_a, sd_b, cov_ab]

    ax.plot(xn, yn, 'r-', label='ODR; '+fit_eq.format(a, b, res.res_var),
            linewidth=2)
    #ax.fill_between(xn, fit_up, fit_dw, alpha=.25, facecolor="r")  # ,label="1-sigma interval")
    ax.legend()
stats.to_csv(opj(idir,'stats_SPM_Rrs.csv'), float_format='%.2f',index=False)
plt.tight_layout(rect=[0.0, 0.0, 0.99, 0.985])
#fig.savefig(opj(figdir,param+'_vs_Rrs.png'), dpi=200)
plt.show()

plt.close()
stats= stats.set_index(['band','algo'])


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

        err = model([stat.b, stat.a], d) - y.values.ravel()

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
fig.savefig(opj(figdir,'histogram_SPM-err-retrieval.png'), dpi=200)
plt.show()

#------------------------------------------------
# compare with uncertainty
#------------------------------------------------
def sig_param(Rrs,sig_a,sig_b,cov_ab):
    return ((Rrs*sig_a)**2 + sig_b**2 + 2*Rrs*cov_ab)**0.5

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

        est = model([stat.b, stat.a], d)
        sig_est = sig_param(d,stat.sig_a,stat.sig_b,stat.cov_ab)
        slope, intercept, r_value, p_value, std_err = scistats.linregress(y, est)
        ax.plot([-10, 100], intercept + slope * np.array([-10, 100]), color=c, lw=1.6)
        ax.errorbar(y, est, xerr=y_sd, yerr=sig_est, fmt='o',
                    color=c, ecolor=c, label=method+'; y={:.1f}x+{:.1f}; $R^2$={:.3f}'.format(slope, intercept,r_value**2),alpha=0.8)



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
fig.savefig(opj(figdir,'compar_SPM-retrieval.png'), dpi=200)
plt.show()

#------------------------------------------------
# simple correlation / linear regression spm/doc
#------------------------------------------------
low = (month >= 9) &(month <= 12)
high = ~low

model,N=linear,2
x=spm
x_sd=spm_sd
y=doc
y_sd=doc_sd
fig, ax = plt.subplots()

ax.errorbar(x[low], y[low], xerr=x_sd[low], yerr=y_sd[low], fmt='o', color='red', ecolor='black', alpha=0.8,label='low')
ax.errorbar(x[high], y[high], xerr=x_sd[high], yerr=y_sd[high], fmt='o', color='blue', ecolor='black', alpha=0.8,label='high')
ax.legend(loc=3)
ax.set_xlabel(r'SSSC (mg/L)')
ax.set_ylabel(r'DOC (mg/L)')

testdata = odr.RealData(x, y, sx=x_sd, sy=y_sd)
_odr = odr.ODR(testdata, odr.Model(model), beta0=[1.] * N)
_odr.set_job(fit_type=2)
res = _odr.run()
res.pprint()
xn = np.linspace(0, np.max(x) * 1.25, 50)
yn = model(res.beta, xn)
fit_up, fit_dw = confidence_interval(xn, model, res)
a, b = res.beta[1], res.beta[0]
ax.plot(xn, yn, 'r-',linewidth=2)
ax.fill_between(xn, fit_up, fit_dw, alpha=.25, facecolor="r")
ax.annotate('y = {:.1f}x + {:.1f}\n$r$ = {:.3f}\nNb = {:d}'.format(
    a, b, np.corrcoef(x, y)[0,1],len(x)),xy=(0.55,0.75), xycoords='axes fraction',)

plt.tight_layout(rect=[0.0, 0.0, 0.99, 0.985])
fig.savefig(opj(figdir,'DOC_vs_SPM.png'), dpi=200)

plt.show()


#------------------------------------------------
# multiparametric fit on wavelength
# SPM = f(Rrs(wl1),...,Rrs(wlN))
#------------------------------------------------
#Rrs := np.array(dims=(Nobs,Nwl))

def multilinear(B, x):
    '''Linear function y = m*x + b'''

    Ndim = x.shape[0]

    fit=B[0]
    for i in range(Ndim):
        fit+= B[i+1] * x[i]
    return fit

Nobs = data.shape[0]
Nwl = 6
Rrs = np.zeros([Nwl,Nobs])
Rrs = data.loc[:, ['Rrs_sim(B3)','Rrs_sim(B4)','Rrs_sim(B5)','Rrs_sim(B6)','Rrs_sim(B7)','Rrs_sim(B8a)']].T
y = data.loc[:, ['TSS (mg L)']]

fitdata = odr.RealData(Rrs, np.array(y).ravel()) #, sx=x_sd, sy=y_sd)
_odr = odr.ODR(fitdata, odr.models.multilinear)
_odr = odr.ODR(fitdata,odr.Model(multilinear), beta0=[1.] * (Nwl+1))
# fit with ordinary least square (OLS, fit_type=2 )
_odr.set_job(fit_type=0)
res = _odr.run()

res.pprint()

for i,SPM in enumerate(y.values):
    print(SPM,multilinear(res.beta,Rrs.iloc[:,i].values))



#------------------------------------------------
# multiparametric fit (and 3D plotting)
#------------------------------------------------

from matplotlib import cm
import scipy.linalg


X,Y = np.meshgrid(np.arange(0, 13, 0.5), np.arange(0,30, 0.5))
XX = X.flatten()
YY = Y.flatten()

# plot Rrs = f(SPM)
nrows=2
ncols=3

fig = plt.figure(figsize=(20, 12))
for i, ax in enumerate(axs.ravel()):
    ax = fig.add_subplot(nrows, ncols, i+1, projection='3d')
    print(i)
    Rrs = data.iloc[:, 2 * i + 2]
    print(Rrs.name)
    Rrs_sd = data.iloc[:, 2 * i + 3]

    data_ = np.array([doc, spm, Rrs]).T

    order = 2    # 1: linear, 2: quadratic
    if order == 1:
        # best-fit linear plane
        A = np.c_[data_[:,0], data_[:,1], np.ones(data_.shape[0])]
        C,_,_,_ = scipy.linalg.lstsq(A, data_[:,2])    # coefficients

        # evaluate it on grid
        Z = C[0]*X + C[1]*Y + C[2]

        # or expressed using matrix/vector product
        #Z = np.dot(np.c_[XX, YY, np.ones(XX.shape)], C).reshape(X.shape)

    elif order == 2:
        # best-fit quadratic curve
        A = np.c_[np.ones(data_.shape[0]), data_[:,:2], np.prod(data_[:,:2], axis=1), data_[:,:2]**2]
        C,_,_,_ = scipy.linalg.lstsq(A, data_[:,2])

        # evaluate it on a grid
        Z = np.dot(np.c_[np.ones(XX.shape), XX, YY, XX*YY, XX**2, YY**2], C).reshape(X.shape)


    ax.plot_surface(X, Y, Z, rstride=1, cstride=1, alpha=0.52, cmap=cm.coolwarm)
    ax.scatter(doc, spm, Rrs, c='r', s=50)
    plt.xlabel('DOC')
    plt.ylabel('SPM')
    ax.set_zlabel('Rrs')

    plt.show()


# functions

# in: [meshgrid], [meshgrid], [np list of coordinate pair np lists. ex: [[x1,y1,z1], [x2,y2,z2], etc.] ], [degree]
# out: [Z]
def curve(X, Y, coord, n):
    XX = X.flatten()
    YY = Y.flatten()

    # best-fit curve
    A = XYchooseN(coord[:,0], coord[:,1], n)
    C,_,_,_ = scipy.linalg.lstsq(A, coord[:,2])
    print("Calculating C in curve() done")
    # evaluate it on a grid
    Z = np.dot(XYchooseN(XX, YY, n), C).reshape(X.shape)
    return Z

# in: [array], [array], [int]
# out: sum from k=0 to k=n of n choose k for x^n-k * y^k (coefficients ignored)
def XYchooseN(x,y,n):
    XYchooseN = []
    n = n+1
    for j in range(len(x)):
        I = x[j]
        J = y[j]
        matrix = []
        Is = []
        Js = []
        for i in range(0,n):
            Is.append(I**i)
            Js.append(J**i)
            matrix.append(np.concatenate((np.ones(n-i),np.zeros(i))))
        Is = np.array(Is)
        Js = np.array(Js)[np.newaxis]
        IsJs0s = matrix * Is * Js.T
        IsJs = []
        for i in range(0,n):
            IsJs = np.concatenate((IsJs,IsJs0s[i,:n-i]))
        XYchooseN.append(IsJs)
    return np.array(XYchooseN)

X,Y = np.meshgrid(x_data,y_data) # surface meshgrid

Z = curve(X, Y, coord, 3)

# todo: cmap by avg from (x1,y1) to (x2,y2) of |Z height - scatter height|

ax.plot_surface(X, Y, Z, rstride=1, cstride=1, alpha=0.5, color='r')




#
#
#
#
# # try band ratio
#     # plot Rrs = f(SPM)
#     fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(30, 12))
#     for i, ax in enumerate(axs.ravel()):
#         print(i)
#         x = data.iloc[:, 2 * i + 2] / data.iloc[:, 2 * 1 + 2]
#
#         print(x.name)
#         x_sd = x * (data.iloc[:, 2 * i + 3] / data.iloc[:, 2 * i + 2] + data.iloc[:, 2 * 1 + 3] / data.iloc[:,
#                                                                                                   2 * 1 + 2])
#
#         ax.errorbar(x, spm, xerr=x_sd, yerr=spm_sd, fmt='o', color='grey', ecolor='black', alpha=0.8)
#         ax.set_title(x.name)
#         ax.set_xlabel(r'$R_{rs}\ (sr^{-1})$')
#         ax.set_ylabel('SPM (mg/L)')
#
#         # fit with ordinary least square (OLS, fit_type=2 )
#         testdata = odr.RealData(x, spm, sx=x_sd, sy=spm_sd)
#         _odr = odr.ODR(testdata, odr.Model(model), beta0=[1.] * N)
#         _odr.set_job(fit_type=2)
#         res = _odr.run()
#         res.pprint()
#         xn = np.linspace(0, np.max(x) * 1.25, 50)
#         yn = model(res.beta, xn)
#         fit_up, fit_dw = confidence_interval(xn, model, res)
#         a, b = res.beta[1], res.beta[0]
#         ax.plot(xn, yn, 'b--', label='OLS; y = {:.1f}x + {:.1f}'.format(a, b), linewidth=2)
#         ax.fill_between(xn, fit_up, fit_dw, alpha=.25, facecolor="b")  # ,label="1-sigma interval")
#
#         # fit with Orthogonal distance regression (ODR, fit_type = 0)
#         _odr.set_job(fit_type=0)
#         res = _odr.run()
#         res.pprint()
#         xn = np.linspace(0, np.max(x) * 1.25, 50)
#         yn = model(res.beta, xn)
#         fit_up, fit_dw = confidence_interval(xn, model, res)
#         a, b = res.beta[1], res.beta[0]
#         ax.plot(xn, yn, 'r-', label='ODR; y = {:.1f}x + {:.1f}'.format(a, b), linewidth=2)
#         ax.fill_between(xn, fit_up, fit_dw, alpha=.25, facecolor="r")  # ,label="1-sigma interval")
#         ax.legend()
#
#
# # test ratio
# x = data.iloc[:, 2 * 5 + 2] / data.iloc[:, 2 * 1 + 2]
