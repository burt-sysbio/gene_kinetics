# read in all files and combine into one large df
import pandas as pd
import os
month = "feb2021"
datadir = "../../data/data_rtm/"
datanames = os.listdir(datadir)

# get filenames for data and corresponding summary files containing fit categories
datanames = [name for name in datanames if ".csv" in name]
studynames = [n[:-4] for n in datanames]

data_list = [pd.read_csv(datadir + n) for n in datanames]

for df,n in zip(data_list, studynames):
    df["ID"] = n

data = pd.concat(data_list)

# need to convert time units
# crawford --> days
# powrie --> days
data["time_unit"] = "hour"
data.loc[data["ID"].str.contains("Craw", regex = False), "time_unit"] = "day"
data.loc[data["ID"].str.contains("Powrie", regex = False), "time_unit"] = "day"
data.loc[data["time_unit"] == "day", "time"] = data.loc[data["time_unit"] == "day", "time"] * 24

# remove isoforms .1 _2 etc from gene names
data[['symbol','iso1', "iso2"]] = data['gene'].str.split(r'_|\.', expand=True)
data.to_csv("../../data/data_summary/data_rtm_combined.csv", index=False)
