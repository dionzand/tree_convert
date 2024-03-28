import json

import pandas as pd

isogg_snp = pd.read_csv("isogg_snps.csv")
unique_subgroups = sorted(list(isogg_snp["Subgroup Name"].unique()))

subsets = []

for subgroup in unique_subgroups:
    for i in range(1, len(subgroup) + 1):
        if subgroup[:i][-1] != "-":
            subsets.append(subgroup[:i])

subsets = sorted(list(set(subsets)))

isogg_tree = {}

for i in subsets:
    if "-" in i or len(i) == 1:
        continue

    if i[:-1] not in isogg_tree:
        isogg_tree[i[:-1]] = [i]
    else:
        isogg_tree[i[:-1]].append(i)

json.dump(isogg_tree, open("isogg_tree2.json", "w"))
