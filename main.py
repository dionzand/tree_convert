import json
import networkx as nx
from collections import Counter
import streamlit as st

isogg_tree_json = json.load(open("isogg_tree.json"))
yfull_tree_json = json.load(open("yfull_tree.json"))

isogg_tree = nx.from_dict_of_lists(isogg_tree_json)
yfull_tree = nx.from_dict_of_lists(yfull_tree_json)

isogg_tree = nx.DiGraph(nx.dfs_tree(isogg_tree, source="ROOT (Y-Chromosome 'Adam')"))
yfull_tree = nx.DiGraph(nx.dfs_tree(yfull_tree, source="ROOT (Y-Chromosome 'Adam')"))

isogg_hg_to_snp_dict = json.load(open("isogg_hg_to_snp_dict.json"))
yfull_hg_to_snp_dict = json.load(open("yfull_hg_to_snp_dict.json"))

isogg_snp_to_hg_dict = json.load(open("isogg_snp_to_hg_dict.json"))
yfull_snp_to_hg_dict = json.load(open("yfull_snp_to_hg_dict.json"))


def get_path_to_root(tree, node):
    path = nx.shortest_path(tree, source="ROOT (Y-Chromosome 'Adam')", target=node)
    return path


def get_identifying_snps_for_path(path, hg_to_snp_dict):
    snps = {}
    for hg in path:
        if hg in hg_to_snp_dict:
            snps[hg] = hg_to_snp_dict[hg]
    return snps


col1, col2 = st.columns(2)

with col1:
    isogg = st.text_input("Enter ISOGG haplogroup")
with col2:
    yfull = st.text_input("Enter YFull haplogroup")

run_button = st.button("Search")

if run_button:
    if "*" in isogg:
        isogg = isogg.split("*")[0]
    if "*" in yfull:
        yfull = yfull.split("*")[0]

    if isogg not in isogg_hg_to_snp_dict:
        st.error(f"ISOGG haplogroup {isogg} not found")
    if yfull not in yfull_hg_to_snp_dict:
        st.error(f"YFull haplogroup {yfull} not found")

    if isogg not in isogg_hg_to_snp_dict or yfull not in yfull_hg_to_snp_dict:
        st.stop()

    else:
        isogg_path = get_path_to_root(isogg_tree, isogg)
        yfull_path = get_path_to_root(yfull_tree, yfull)

        isogg_identifying_snps = [i for i in isogg_hg_to_snp_dict.get(isogg)]  # if i in yfull_snp_to_hg_dict
        yfull_identifying_snps = [i for i in yfull_hg_to_snp_dict.get(yfull)]  # if i in isogg_snp_to_hg_dict

        isogg_missing_snps = [i for i in isogg_identifying_snps if i not in yfull_snp_to_hg_dict]
        yfull_missing_snps = [i for i in yfull_identifying_snps if i not in isogg_snp_to_hg_dict]

        if len(isogg_missing_snps) > 0:
            st.error(f"These ISOGG SNPs are not present in YFull: {', '.join(isogg_missing_snps)}")
        if len(yfull_missing_snps) > 0:
            st.error(f"These YFull SNPs are not present in ISOGG: {', '.join(yfull_missing_snps)}")

        isogg_present_snps = [i for i in isogg_identifying_snps if i in yfull_snp_to_hg_dict]
        yfull_present_snps = [i for i in yfull_identifying_snps if i in isogg_snp_to_hg_dict]

        yfull_hg_isogg_snps = Counter([i for snp in isogg_present_snps for i in yfull_snp_to_hg_dict[snp]])
        isogg_hg_yfull_snps = Counter([i for snp in yfull_present_snps for i in isogg_snp_to_hg_dict[snp]])

        st.header("ISOGG resolution")
        st.write(f"**There are {len(isogg_hg_yfull_snps)} ISOGG haplogroups that correspond to these YFull SNPs:**")
        for hg in isogg_hg_yfull_snps:
            st.write(f"Comparing Yfull based ISOGG haplogroup {hg} ({int(isogg_hg_yfull_snps[hg]/sum(isogg_hg_yfull_snps.values())*100)}%) to real ISOGG haplogroup {isogg}:")
            path = get_path_to_root(isogg_tree, hg)
            if path == isogg_path:
                st.warning("Resolutions are the same")
            elif len(path) < len(isogg_path):
                if path == isogg_path[:len(path)]:
                    st.success(f"ISOGG resolution is {len(isogg_path) - len(path)} levels higher")
                else:
                    lowest_common_ancestor = nx.lowest_common_ancestor(isogg_tree, path[-1], isogg_path[-1])
                    st.error(f"Paths do not match. Lowest common ancestor is {lowest_common_ancestor}.")
            elif len(path) > len(isogg_path):
                if isogg_path == path[:len(isogg_path)]:
                    st.success(f"YFull resolution is {len(path) - len(isogg_path)} levels higher")
                else:
                    lowest_common_ancestor = nx.lowest_common_ancestor(isogg_tree, path[-1], isogg_path[-1])
                    st.error(f"Paths do not match. Lowest common ancestor is {lowest_common_ancestor}.")

        st.header("YFull resolution")
        st.write(f"**There are {len(yfull_hg_isogg_snps)} YFull haplogroups that correspond to these ISOGG SNPs:**")
        for hg in yfull_hg_isogg_snps:
            st.write(f"Comparing ISOGG based YFull haplogroup {hg} ({int(yfull_hg_isogg_snps[hg]/sum(yfull_hg_isogg_snps.values())*100)}%) to real YFull haplogroup {yfull}:")
            path = get_path_to_root(yfull_tree, hg)
            if path == yfull_path:
                st.warning("Resolutions are the same")
            elif len(path) < len(yfull_path):
                if path == yfull_path[:len(path)]:
                    st.success(f"YFull resolution is {len(yfull_path) - len(path)} levels higher")
                else:
                    lowest_common_ancestor = nx.lowest_common_ancestor(yfull_tree, path[-1], yfull_path[-1])
                    st.error(f"Paths do not match. Lowest common ancestor is {lowest_common_ancestor}.")
            elif len(path) > len(yfull_path):
                if yfull_path == path[:len(yfull_path)]:
                    st.success(f"ISOGG resolution is {len(path) - len(yfull_path)} levels higher")
                else:
                    lowest_common_ancestor = nx.lowest_common_ancestor(yfull_tree, path[-1], yfull_path[-1])
                    st.error(f"Paths do not match. Lowest common ancestor is {lowest_common_ancestor}.")

        st.header("Paths")
        st.write(f"Yfull path is:")
        st.info(' - '.join(yfull_path))
        st.write(f"ISOGG path is:")
        st.info(' - '.join(isogg_path))

        st.header("Identifying SNPs")
        st.write(f"Identifying SNPs for Yfull {yfull} are:")
        st.info(', '.join(yfull_identifying_snps))
        st.write(f"Identifying SNPs for ISOGG {isogg} are:")
        st.info(', '.join(isogg_identifying_snps))

        st.header("Haplogroups")
        st.write(f"ISOGG haplogroups for these Yfull SNPs are: {isogg_hg_yfull_snps}")
        st.write(f"YFull haplogroups for these ISOGG SNPs are: {yfull_hg_isogg_snps}")

        st.header("Info")
        st.info("The resolution of the real ISOGG haplogroup is compared to the resolution of the YFull haplogroup based on ISOGGs identifying SNPs. "
                "The percentage of each haplogroup is based on the number of identifying SNPs that support this haplogroup. "
                "Asterisks (*) are currently ignored for comparison.")
