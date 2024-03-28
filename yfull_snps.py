import json
import pandas as pd

yfull_snps = pd.read_csv("yfull_snps.csv")

print(yfull_snps.head())

yfull_snp_to_hg_dict = {}
yfull_hg_to_snp_dict = {}

# convert to dict where key is Name and value is Subgroup Name
for index, row in yfull_snps.iterrows():
    snp = row["name"].replace("^", "").replace("^^", "")
    hg = row["name.1"]
    if hg not in yfull_hg_to_snp_dict:
        yfull_hg_to_snp_dict[hg] = [snp]
    else:
        yfull_hg_to_snp_dict[hg].append(snp)

    if snp not in yfull_snp_to_hg_dict:
        yfull_snp_to_hg_dict[snp] = [hg]
    else:
        yfull_snp_to_hg_dict[snp].append(hg)

# save to json
json.dump(yfull_snp_to_hg_dict, open("yfull_snp_to_hg_dict.json", "w"))
json.dump(yfull_hg_to_snp_dict, open("yfull_hg_to_snp_dict.json", "w"))
