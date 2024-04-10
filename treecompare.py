import argparse
import json

import networkx as nx
from collections import Counter

import pandas as pd

isogg_tree_json = json.load(open("data/isogg_tree.json"))
yfull_tree_json = json.load(open("data/yfull_tree.json"))

isogg_tree = nx.from_dict_of_lists(isogg_tree_json)
yfull_tree = nx.from_dict_of_lists(yfull_tree_json)

isogg_tree = nx.DiGraph(nx.dfs_tree(isogg_tree, source="ROOT (Y-Chromosome 'Adam')"))
yfull_tree = nx.DiGraph(nx.dfs_tree(yfull_tree, source="ROOT (Y-Chromosome 'Adam')"))

isogg_hg_to_snp_dict = json.load(open("data/isogg_hg_to_snp_dict.json"))
yfull_hg_to_snp_dict = json.load(open("data/yfull_hg_to_snp_dict.json"))

isogg_snp_to_hg_dict = json.load(open("data/isogg_snp_to_hg_dict.json"))
yfull_snp_to_hg_dict = json.load(open("data/yfull_snp_to_hg_dict.json"))


def get_path_to_root(tree, node):
    path = nx.shortest_path(tree, source="ROOT (Y-Chromosome 'Adam')", target=node)
    return path


def get_identifying_snps_for_path(path, hg_to_snp_dict):
    snps = {}
    for hg in path:
        if hg in hg_to_snp_dict:
            snps[hg] = hg_to_snp_dict[hg]
    return snps


def main_cli(sample, isogg, yfull):
    log = []

    # TODO: Add exlusion of * from haplogroups
    if "*" in isogg:
        isogg = isogg.split("*")[0]
    if "*" in yfull:
        yfull = yfull.split("*")[0]

    if isogg not in isogg_hg_to_snp_dict:
        log.append(f"No known SNPs for ISOGG haplogroup {isogg}")
    if yfull not in yfull_hg_to_snp_dict:
        log.append(f"No known SNPs for YFull haplogroup {yfull}")

    if isogg not in isogg_hg_to_snp_dict or yfull not in yfull_hg_to_snp_dict:
        return {
            "sample": sample,
            "isogg": isogg,
            "yfull": yfull,
            "log": log,
        }

    else:
        isogg_path = get_path_to_root(isogg_tree, isogg)
        yfull_path = get_path_to_root(yfull_tree, yfull)

        isogg_identifying_snps = [i for i in isogg_hg_to_snp_dict.get(isogg)]
        yfull_identifying_snps = [i for i in yfull_hg_to_snp_dict.get(yfull)]

        isogg_missing_snps = [i for i in isogg_identifying_snps if i not in yfull_snp_to_hg_dict]
        yfull_missing_snps = [i for i in yfull_identifying_snps if i not in isogg_snp_to_hg_dict]

        isogg_present_snps = [i for i in isogg_identifying_snps if i in yfull_snp_to_hg_dict]
        yfull_present_snps = [i for i in yfull_identifying_snps if i in isogg_snp_to_hg_dict]

        yfull_hg_isogg_snps = Counter([i for snp in isogg_present_snps for i in yfull_snp_to_hg_dict[snp]])
        isogg_hg_yfull_snps = Counter([i for snp in yfull_present_snps for i in isogg_snp_to_hg_dict[snp]])

        if len(isogg_hg_yfull_snps) == 0:
            log.append("No ISOGG haplogroups found for these YFull SNPs")
            highest_isogg_hg = None
            highest_isogg_hg_ratio = None
            isogg_resolution = None
        else:
            highest_isogg_hg = max(isogg_hg_yfull_snps, key=isogg_hg_yfull_snps.get)
            if list(isogg_hg_yfull_snps.values()).count(isogg_hg_yfull_snps[highest_isogg_hg]) > 1:
                log.append(f"Warning: Multiple ISOGG haplogroups share the same highest value. Might be ambiguous.")
                # If multiple haplogroups share the same highest value, choose the one with the longest path to the root
                highest_isogg_hg = max([hg for hg in isogg_hg_yfull_snps if isogg_hg_yfull_snps[hg] == isogg_hg_yfull_snps[highest_isogg_hg]],
                                       key=lambda hg: len(get_path_to_root(isogg_tree, hg)))
            highest_isogg_hg_ratio = isogg_hg_yfull_snps[highest_isogg_hg] / sum(isogg_hg_yfull_snps.values())
            highest_isogg_path = get_path_to_root(isogg_tree, highest_isogg_hg)
            if highest_isogg_path == isogg_path:
                isogg_resolution = "Resolutions are the same"
            elif len(highest_isogg_path) < len(isogg_path):
                if highest_isogg_path == isogg_path[:len(highest_isogg_path)]:
                    isogg_resolution = f"ISOGG resolution is {len(isogg_path) - len(highest_isogg_path)} levels higher"
                else:
                    lowest_common_ancestor = nx.lowest_common_ancestor(isogg_tree, highest_isogg_path[-1], isogg_path[-1])
                    isogg_resolution = f"Paths do not match. Lowest common ancestor is {lowest_common_ancestor}."
            elif len(highest_isogg_path) > len(isogg_path):
                if isogg_path == highest_isogg_path[:len(isogg_path)]:
                    isogg_resolution = f"YFull resolution is {len(highest_isogg_path) - len(isogg_path)} levels higher"
                else:
                    lowest_common_ancestor = nx.lowest_common_ancestor(isogg_tree, highest_isogg_path[-1], isogg_path[-1])
                    isogg_resolution = f"Paths do not match. Lowest common ancestor is {lowest_common_ancestor}."
            else:
                lowest_common_ancestor = nx.lowest_common_ancestor(isogg_tree, highest_isogg_path[-1], isogg_path[-1])
                isogg_resolution = f"Paths do not match. Lowest common ancestor is {lowest_common_ancestor}."

        if len(yfull_hg_isogg_snps) == 0:
            log.append("No YFull haplogroups found for these ISOGG SNPs")
            highest_yfull_hg = None
            highest_yfull_hg_ratio = None
            yfull_resolution = None
        else:
            highest_yfull_hg = max(yfull_hg_isogg_snps, key=yfull_hg_isogg_snps.get)
            if list(yfull_hg_isogg_snps.values()).count(yfull_hg_isogg_snps[highest_yfull_hg]) > 1:
                log.append(f"Warning: Multiple YFull haplogroups share the same highest value. Might be ambiguous.")
                highest_yfull_hg = max([hg for hg in yfull_hg_isogg_snps if yfull_hg_isogg_snps[hg] == yfull_hg_isogg_snps[highest_yfull_hg]],
                                       key=lambda hg: len(get_path_to_root(yfull_tree, hg)))
            highest_yfull_hg_ratio = yfull_hg_isogg_snps[highest_yfull_hg] / sum(yfull_hg_isogg_snps.values())
            highest_yfull_path = get_path_to_root(yfull_tree, highest_yfull_hg)

            if highest_yfull_path == yfull_path:
                yfull_resolution = "Resolutions are the same"
            elif len(highest_yfull_path) < len(yfull_path):
                if highest_yfull_path == yfull_path[:len(highest_yfull_path)]:
                    yfull_resolution = f"YFull resolution is {len(yfull_path) - len(highest_yfull_path)} levels higher"
                else:
                    lowest_common_ancestor = nx.lowest_common_ancestor(yfull_tree, highest_yfull_path[-1], yfull_path[-1])
                    yfull_resolution = f"Paths do not match. Lowest common ancestor is {lowest_common_ancestor}."
            elif len(highest_yfull_path) > len(yfull_path):
                if yfull_path == highest_yfull_path[:len(yfull_path)]:
                    yfull_resolution = f"ISOGG resolution is {len(highest_yfull_path) - len(yfull_path)} levels higher"
                else:
                    lowest_common_ancestor = nx.lowest_common_ancestor(yfull_tree, highest_yfull_path[-1], yfull_path[-1])
                    yfull_resolution = f"Paths do not match. Lowest common ancestor is {lowest_common_ancestor}."
            else:
                lowest_common_ancestor = nx.lowest_common_ancestor(yfull_tree, highest_yfull_path[-1], yfull_path[-1])
                yfull_resolution = f"Paths do not match. Lowest common ancestor is {lowest_common_ancestor}."

        return {
            "sample": sample,
            "isogg": isogg,
            "yfull": yfull,
            "isogg_path": isogg_path,
            "yfull_path": yfull_path,
            "isogg_identifying_snps": isogg_identifying_snps,
            "yfull_identifying_snps": yfull_identifying_snps,
            "isogg_missing_snps": isogg_missing_snps,
            "yfull_missing_snps": yfull_missing_snps,
            "isogg_hg_yfull_snps": isogg_hg_yfull_snps,
            "yfull_hg_isogg_snps": yfull_hg_isogg_snps,
            "highest_isogg_hg": highest_isogg_hg,
            "highest_isogg_hg_ratio": highest_isogg_hg_ratio,
            "highest_yfull_hg": highest_yfull_hg,
            "highest_yfull_hg_ratio": highest_yfull_hg_ratio,
            "isogg_resolution": isogg_resolution,
            "yfull_resolution": yfull_resolution,
            "log": log,
        }


if __name__ == '__main__':
    argparser = argparse.ArgumentParser()
    argparser.add_argument("--input",
                           help="Input file with haplogroups to compare (in CSV format, "
                                "with three columns: Sample, ISOGG_hg, YFull_hg)",
                           required=True)
    argparser.add_argument("--output",
                           help="Output file for results (in CSV format)",
                           default="results.csv")
    args = argparser.parse_args()

    input_csv = pd.read_csv(args.input, names=["Sample", "ISOGG_hg", "YFull_hg"], header=0)

    results_list = []

    for _, row in input_csv.iterrows():
        result = main_cli(row["Sample"], row["ISOGG_hg"], row["YFull_hg"])
        results_list.append(result)

    results_df = pd.DataFrame(results_list)
    results_df.to_csv(args.output, index=False)
