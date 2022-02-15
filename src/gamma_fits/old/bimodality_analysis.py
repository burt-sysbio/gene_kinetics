
import pandas as pd
import numpy as np
from statsmodels.distributions.empirical_distribution import ECDF
from scipy.stats import gamma, uniform


def simulate_dip(data, pdftype = "uniform", nsim = 1000, signif = 99):
    """
    data should be df specific for one gene and one celltype and one study
    signif is significance level, use 95 for pval 0.05 and 99 for pval 0.01
    """
    # number of unique timepoints
    nobs = data["nobs"]
    alpha = data["alpha"]
    beta = data["beta"]

    # initiate a pdf with given parameters
    if pdftype == "gamma":
        func = gamma(a=alpha, scale = 1/beta)
    elif pdftype == "uniform":
        func = uniform(scale = data["time"])

    def get_dip(nobs, func):
        # draw samples for comparison dip
        sample = func.rvs(size = nobs)

        # compute empirical distribtion function for samples
        ecdf_vals = ECDF(sample)

        # get the true y values from dist based on the x values from ecdf
        y = func.cdf(ecdf_vals.x)

        # get dip
        dip = max(np.abs(y - ecdf_vals.y))

        return dip

    dip_list = [get_dip(nobs, func) for _ in range(nsim)]
    dip_list = np.asarray(dip_list)

    perc = np.percentile(dip_list, 99)
    return dip_list, perc


def get_dip(data):
    """
    purpose: plot gene expression time course from data with corresponding fits
    data: df from data_rtm processed data
    fit : df from gamma fits generated from data
    """
    # get fit vals but only until 1 is reached

    x_idx = data.avg_norm_rtm2.values
    x_idx = np.where(x_idx == 1)[0][0]
    data = data.iloc[:(x_idx+1), :]

    alpha = data["alpha"].iloc[0]
    beta = data["beta"].iloc[0]

    # get gamma cdf vals
    gamma_ydata = gamma.cdf(data["time"], a = alpha, scale = 1 / beta)

    dip = max(np.abs(data["avg_norm_rtm2"] - gamma_ydata))

    return dip


# read in all files
month = "feb2021"
fitdir = "../../output/gamma_fits/" + month + "/"
datadir = "../../data/data_rtm/"


studies = [
    "Crawford_rtm_CD4.arm_lcmv",
    "Crawford_rtm_CD4.cl13_lcmv",
    "Crawford_rtm_CD8.arm_lcmv",
    "Crawford_rtm_CD8.cl13_lcmv",
    "Powrie_rtm_innate_colitis",
    "Proserpio_rtm_Th2_parasite",
    "Peine_rtm_Th0_invitro",
    "Peine_rtm_Th1_invitro",
    "Peine_rtm_Th2_invitro",
    "Peine_rtm_ThMix_invitro",
    "Nir_rtm_Th0_invitro",
    "Nir_rtm_Th17_invitro"
]
#study = studies[1]
#study = "Proserpio_rtm_Th2_parasite.csv"

bimodal_studies = []
for study in studies:
    print("running diptest for " + study)
    fFit = "fit_res_" + study + ".csv"
    fdata = study + ".csv"
    fGenes = "fit_summary_" + study + ".csv"

    fit_res = pd.read_csv(fitdir + fFit)
    data = pd.read_csv(datadir + fdata)
    df_genes = pd.read_csv(fitdir + fGenes)

    # get all genes in other category and run bimodality test on them
    df_other = df_genes.loc[df_genes["best_fit"] == "other", "gene"]
    fit_other = fit_res.loc[fit_res["gene"].isin(df_other)]

    # only focus on gamma fit for dip test...?
    fit_other = fit_other.loc[fit_other["model"] == "gamma"]

    # get number of observations from data, needed to sample monte carlo dips from uniform
    df_obs = data.groupby('gene')['time'].nunique()
    df_obs = df_obs.reset_index()
    df_obs.columns = ["gene", "nobs"]

    # merge df fit with df that has number of data observations
    df_diptest = pd.merge(fit_other, df_obs, how = "left", on = "gene")

    # add the last time point to the data to get scale for uniform dist
    df_lasttime = data.groupby(["gene"])["time"].max().reset_index()

    # combine
    df_diptest = pd.merge(df_diptest, df_lasttime, how = "left", on = "gene")

    # add alpha and beta parameters to the data df
    data_other = data.loc[data["gene"].isin(fit_other["gene"])]

    data_other = pd.merge(data_other, fit_other[["gene", "alpha", "beta"]], on = "gene", how = "left")

    # sample some genes
    genes = df_diptest["gene"].drop_duplicates()

    # for each gene if it is bimodal, attach
    bimodal_genes = []
    for gene in genes:

        gene_diptest = df_diptest.loc[df_diptest["gene"] == gene, :]
        diplist, perc = simulate_dip(gene_diptest)

        gene_data = data_other.loc[data_other["gene"] == gene, :]
        dip_data = get_dip(gene_data)

        if dip_data > perc:
            bimodal_genes.append(gene)

    bimodal_studies.append(bimodal_genes)

# store output
df_bimodal = pd.DataFrame(bimodal_studies).T
df_bimodal.columns = studies
df_bimodal.to_csv("../../output/bimodal_genes/bimodal_genes_thres001.csv", index = False)
