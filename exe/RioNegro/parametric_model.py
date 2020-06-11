import numpy as np
import pandas as pd

import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rcParams.update({'font.size': 18})
import statsmodels.api as sm
from sklearn.linear_model import LinearRegression
from statsmodels.sandbox.regression.predstd import wls_prediction_std
from scipy import odr as odr

file = '/DATA/OBS2CO/data/rogerio/data/Simulation_Rrs_OSOAA_TSS.xlsx'
data = pd.read_excel(file)
data.sort_values(by="TSS (mg L)", inplace=True)
data = data.set_index(['Station', 'Date'])
spm = data.iloc[:, 0]  # .values.reshape(-1, 1)
spm_sd = data.iloc[:, 1]


def linear(B, x):
    '''Linear function y = m*x + b'''
    return B[0] + B[1] * x


def exponential(B, x):
    return B[0] * np.exp(B[1] * x)


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


models = [(linear, 2)]  # (poly,3),,(exponential,2)
for model, N in models:

    # plot Rrs = f(SPM)
    fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(25, 12))
    for i, ax in enumerate(axs.ravel()):
        print(i)
        x = data.iloc[:, 2 * i + 2]
        print(x.name)
        x_sd = data.iloc[:, 2 * i + 3]

        ax.errorbar(x, spm, xerr=x_sd, yerr=spm_sd, fmt='o', color='grey', ecolor='black', alpha=0.8)
        ax.set_title(x.name)
        ax.set_xlabel(r'$R_{rs}\ (sr^{-1})$')
        ax.set_ylabel('SPM (mg/L)')

        testdata = odr.RealData(x, spm, sx=x_sd, sy=spm_sd)
        _odr = odr.ODR(testdata, odr.Model(model), beta0=[1.] * N)

        # fit with ordinary least square (OLS, fit_type=2 )
        _odr.set_job(fit_type=2)
        res = _odr.run()
        res.pprint()
        xn = np.linspace(0, np.max(x) * 1.25, 50)
        yn = model(res.beta, xn)
        fit_up, fit_dw = confidence_interval(xn, model, res)
        a, b = res.beta[1], res.beta[0]
        ax.plot(xn, yn, 'b--', label='OLS; y = {:.1f}x + {:.1f}\n res_var = {:.1f}'.format(a, b, res.res_var),
                linewidth=2)
        ax.fill_between(xn, fit_up, fit_dw, alpha=.25, facecolor="b")  # ,label="1-sigma interval")

        # fit with Orthogonal distance regression (ODR, fit_type = 0)
        _odr.set_job(fit_type=0)
        res = _odr.run()
        res.pprint()
        xn = np.linspace(0, np.max(x) * 1.25, 50)
        yn = model(res.beta, xn)
        fit_up, fit_dw = confidence_interval(xn, model, res)
        a, b = res.beta[1], res.beta[0]
        ax.plot(xn, yn, 'r-', label='ODR; y = {:.1f}x + {:.1f}\n res_var = {:.1f}'.format(a, b, res.res_var),
                linewidth=2)
        ax.fill_between(xn, fit_up, fit_dw, alpha=.25, facecolor="r")  # ,label="1-sigma interval")
        ax.legend()
    plt.tight_layout(rect=[0.0, 0.0, 0.99, 0.985])
    fig.savefig('/DATA/ARTICLES/rogerio/SPM_vs_Rrs.png', dpi=200)

    # plot Rrs = f(SPM)
    fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(30, 12))
    for i, ax in enumerate(axs.ravel()):
        print(i)
        x = data.iloc[:, 2 * i + 2] / data.iloc[:, 2 * 1 + 2]

        print(x.name)
        x_sd = x * (data.iloc[:, 2 * i + 3] / data.iloc[:, 2 * i + 2] + data.iloc[:, 2 * 1 + 3] / data.iloc[:,
                                                                                                  2 * 1 + 2])

        ax.errorbar(x, spm, xerr=x_sd, yerr=spm_sd, fmt='o', color='grey', ecolor='black', alpha=0.8)
        ax.set_title(x.name)
        ax.set_xlabel(r'$R_{rs}\ (sr^{-1})$')
        ax.set_ylabel('SPM (mg/L)')

        # fit with ordinary least square (OLS, fit_type=2 )
        testdata = odr.RealData(x, spm, sx=x_sd, sy=spm_sd)
        _odr = odr.ODR(testdata, odr.Model(model), beta0=[1.] * N)
        _odr.set_job(fit_type=2)
        res = _odr.run()
        res.pprint()
        xn = np.linspace(0, np.max(x) * 1.25, 50)
        yn = model(res.beta, xn)
        fit_up, fit_dw = confidence_interval(xn, model, res)
        a, b = res.beta[1], res.beta[0]
        ax.plot(xn, yn, 'b--', label='OLS; y = {:.1f}x + {:.1f}'.format(a, b), linewidth=2)
        ax.fill_between(xn, fit_up, fit_dw, alpha=.25, facecolor="b")  # ,label="1-sigma interval")

        # fit with Orthogonal distance regression (ODR, fit_type = 0)
        _odr.set_job(fit_type=0)
        res = _odr.run()
        res.pprint()
        xn = np.linspace(0, np.max(x) * 1.25, 50)
        yn = model(res.beta, xn)
        fit_up, fit_dw = confidence_interval(xn, model, res)
        a, b = res.beta[1], res.beta[0]
        ax.plot(xn, yn, 'r-', label='ODR; y = {:.1f}x + {:.1f}'.format(a, b), linewidth=2)
        ax.fill_between(xn, fit_up, fit_dw, alpha=.25, facecolor="r")  # ,label="1-sigma interval")
        ax.legend()


# test ratio
x = data.iloc[:, 2 * 5 + 2] / data.iloc[:, 2 * 1 + 2]
