from scipy.special import gammainc
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.special import gammainc
from scipy import stats
from statsmodels.stats.multitest import fdrcorrection


def gamma_cdf(t, alpha, beta):
    dummy = t * beta
    return gammainc(alpha, dummy)


def gamma_cdf1(t, beta):
    return gamma_cdf(t, 1, beta)


def conv_cols(df):
    df.columns = df.columns.droplevel()
    df.columns = [str(col) for col in df.columns]
    df = df.reset_index()
    return (df)


def prep_data(df):
    # check column names because old version used gene name
    if "gene" in df.columns:
        df = df.rename(columns = {"gene" : "gene_name"})

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
        # check that there are no nans in y
        assert (~np.isnan(y).any())
        # check that data is normalized and keep array only until max is reached
        assert np.max(y) == 1.0
        max_idx = np.where(y == 1.0)[0][0]
        y = y[:(max_idx + 1)]

        # only focus on arrays where y has 4 time points
        if len(y) > 3:
            xdata = x[:len(y)]
            if errs.size != 0:
                sigma = errs[i, :]
                # sigma = sigma[~np.isnan(sigma)]
                sigma = sigma[:max_idx]
                sigma_abs = True
            else:
                sigma = None
                sigma_abs = False

            # run fit and append output including gene name
            out = fit_gamma(xdata, y, sigma, sigma_abs, gamma_fun)
            out = [gene] + out
            fit_res.append(out)

    # convert fit res to dataframe
    colnames = ["gene", "alpha", "beta", "rss", "alpha_err", "beta_err"]
    df_fit_res = pd.DataFrame(fit_res, columns= colnames)
    return df_fit_res


def fit_gamma(x, y, sigma, sigma_abs, gamma_fun):
    # dummy output in case fit does not work (e.g. gene kinetics are bimodal)
    out = [np.nan, np.nan, np.nan]
    try:
        # run fit catch runtime exception for bad fit
        fit_val, fit_err = curve_fit(f=gamma_fun, xdata=x, ydata=y, sigma=sigma,
                                     absolute_sigma=sigma_abs)

        # error of fitted params (not used atm)
        err = np.sqrt(np.diag(fit_err))

        # assign fit values depending on gamma_fun (exponential fit only fits beta not alpha)
        if gamma_fun == gamma_cdf:
            alpha_fit, beta_fit = fit_val
            alpha_err, beta_err = err
        else:
            beta_fit = fit_val[0]
            alpha_fit = 1
            alpha_err = np.nan
            beta_err = err[0]

        # compute chi sqr and get residuals
        # need to change pipeline function in R to write sigma as 1 instead of NA
        nexp = gamma_fun(x, *fit_val)
        r = y - nexp
        if sigma is None: sigma = 1
        rss = np.sum((r / sigma) ** 2)

        out = [alpha_fit, beta_fit, rss, alpha_err, beta_err]

    except RuntimeError:
        print("max calls reached no good fit")

    return out


def f_test(rss1, rss2, n_obs, f1=1, f2=2):
    d1 = f2 - f1
    d2 = n_obs - f2
    F = ((rss1 - rss2) / d1) / (rss2 / d2)
    f_crit = 1 - stats.f.cdf(F, d1, d2)
    return f_crit


def run_f_test(df, gamma_1=gamma_cdf1, gamma_2=gamma_cdf):
    fit_res1 = fit_kinetic(df, gamma_1)
    fit_res2 = fit_kinetic(df, gamma_2)

    df_fit_res = pd.merge(fit_res1, fit_res2, on="gene")

    # get genes were any fit did not work well
    genes_bimodal = np.isnan(df_fit_res[["alpha_x", "alpha_y"]]).any(axis=1)
    df_bimodal = df_fit_res.loc[genes_bimodal]
    df_fit_res = df_fit_res.loc[~genes_bimodal]

    # run f test to compare gamma fits
    timepoints = df.time.drop_duplicates().values
    n_timepoints = len(timepoints)
    df_fit_res["pval"] = f_test(df_fit_res["rss_x"].values,
                                df_fit_res["rss_y"].values,
                                n_obs=n_timepoints)

    # get fdr, use only second element which are adjusted pvalues
    df_fit_res["padj"] = fdrcorrection(df_fit_res.pval.values)[1]
    df_fit_res["sig"] = "Non-exponential"
    df_fit_res["sig"][df_fit_res["padj"]>0.05] = "Exponential"

    # add genes for which fit did not work
    df_bimodal["pval"] = None
    df_bimodal["padj"] = None
    df_bimodal["sig"] = "other"

    df_fit_res = pd.concat([df_fit_res, df_bimodal])
    return df_fit_res

#df = pd.read_csv("../../output/rtm_data/peine_rtm_Th0_invitro.csv")
#fit_res = run_f_test(df)