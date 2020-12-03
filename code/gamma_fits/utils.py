from scipy.special import gammainc
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.special import gammainc

def gamma_cdf(t, alpha, beta):
    dummy = t*beta
    return gammainc(alpha,dummy)

def gamma_cdf1(t, beta):
    dummy = t*beta
    return gammainc(1.0,dummy)


def gamma_cdf2(t, beta):
    dummy = t*beta
    return gammainc(2,dummy)


def gamma_cdf3(t, beta):
    dummy = t*beta
    return gammainc(10.0,dummy)


def conv_cols(df):
    df.columns = df.columns.droplevel()
    df.columns = [str(col) for col in df.columns]
    df = df.reset_index()
    return (df)


def prep_data(df):
    df_err = df[["cell_type", "gene_name", "time", "SD"]]
    df_err = df_err.drop_duplicates()

    df_val = df[["cell_type", "gene_name", "time", "avg_norm_rtm2"]]
    df_val = df_val.drop_duplicates()

    # need to make wide again for fit
    df_list = [df_val, df_err]
    df_list = [pd.pivot_table(data=df,
                              index=["cell_type", "gene_name"],
                              columns="time") for df in df_list]
    df_list = [conv_cols(df) for df in df_list]

    # extract data and errors
    timepoints = df.time.drop_duplicates().values
    n_timepoints = len(timepoints)

    genes = df_list[0].gene_name.values
    vals = df_list[0].iloc[:, -n_timepoints:].values
    errs = df_list[1].iloc[:, -n_timepoints:].values

    return genes, vals, errs, timepoints


def fit_kinetic(df, gamma_fun):
    genes, vals, errs, x = prep_data(df)

    fit_res = []
    # fit gamma dist for each gene
    for i, gene in zip(range(len(vals)), genes):
        y = vals[i, :]
        #check that data is normalized and keep array only until max is reached
        assert np.max(y) == 1.0
        max_idx = np.where(y == 1.0)[0][0]
        y = y[:max_idx]
        # kick out nans in y and later in sigma
        # y = y[~np.isnan(y)]
        # only focus on arrays where y has 4 time points
        if len(y) > 3:
            xdata = x[:len(y)]
            if errs.size != 0:
                sigma = errs[i, :]
                # sigma = sigma[~np.isnan(sigma)]
                sigma = sigma[:max_idx]
                absolute_sigma = True
            else:
                sigma = None
                absolute_sigma = False
            alpha_fit, beta_fit, chisqr = fit_gamma(xdata, y, sigma, absolute_sigma, gamma_fun)
            out = [gene, alpha_fit, beta_fit, chisqr]
        # add alpha fit and err
        fit_res.append(out)

    return fit_res


def fit_gamma(x, y, sigma, absolute_sigma, gamma_fun):
    # check that data is normalized and get index of max, check that at least 4 legit time points are there
    # alpha_fit = np.nan
    chisq = np.nan
    out = np.array([0, 0, 0])
    # only do fit procedure if at least 3 vals in y to fit
    try:
        # run fit catch runtime exception for bad fit
        # restrain alpha within 1 and 100
        fit_val, fit_err = curve_fit(f=gamma_fun,xdata=x, ydata=y,sigma=sigma,absolute_sigma=absolute_sigma)
        # according to docs this is the error of covariance of fitted params
        alpha_err = np.sqrt(np.diag(fit_err))

        if gamma_fun == gamma_cdf:
            alpha_fit, beta_fit = fit_val
        else:
            beta_fit = fit_val
            alpha_fit = 1

        # compute chi sqr
        nexp = gamma_fun(x, *fit_val)
        # get residuals
        r = y - nexp
        # only reassign chisq for low error in identified parameter
        if (alpha_err < 10.0).all():
            if sigma is None:
                chisq = np.sum(r ** 2)
            else:
                chisq = np.sum((r / sigma) ** 2)

        out = np.array([alpha_fit, beta_fit, chisq])

    except RuntimeError:
        print("max calls reached no good fit")
        print(x)
        print(y)

    return out


def run_fit(df, data, cell, inf, gamma_list=[gamma_cdf, gamma_cdf1]):

    df_list, fit_res = fit_kinetic(df, gamma_list=gamma_list)
    #df_final = df_list[0]
    #df_final["alpha1"] = err_list[0]
    #df_final["alpha2"] = err_list[1]
    #df_final["alpha10"] = err_list[2]

    #df_final["beta_1"] = beta_list[0]
    #df_final["beta_2"] = beta_list[1]
    #df_final["beta_10"] = beta_list[2]

    # nas were generated during fit if not enough data was available
    #df_final = df_final.dropna()
    # df_final = df_final[["cell_type", "gene_name", "alpha1", "alpha2", "alpha10"]]
    # df_final = df_final.reset_index(drop=True)

    #filename = "output/gamma_fits/gamma_fit_" + data + "_" + cell + "_" + inf + ".csv"
    #df_final.to_csv(filename, index=False)

    return df_final
