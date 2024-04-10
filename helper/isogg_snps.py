import json
import pandas as pd

isogg_snp = pd.read_csv("../data/isogg_snps.csv")

isogg_snp_to_hg_dict = {}
isogg_hg_to_snp_dict = {}

# convert to dict where key is Name and value is Subgroup Name
for index, row in isogg_snp.iterrows():
    snp = row["Name"].replace("^", "").replace("^^", "")
    hg = row["Subgroup Name"]
    if hg not in isogg_hg_to_snp_dict:
        isogg_hg_to_snp_dict[hg] = [snp]
    else:
        isogg_hg_to_snp_dict[hg].append(snp)

    if snp not in isogg_snp_to_hg_dict:
        isogg_snp_to_hg_dict[snp] = [hg]
    else:
        isogg_snp_to_hg_dict[snp].append(hg)

# save to json
json.dump(isogg_snp_to_hg_dict, open("../data/isogg_snp_to_hg_dict.json", "w"))
json.dump(isogg_hg_to_snp_dict, open("../data/isogg_hg_to_snp_dict.json", "w"))
