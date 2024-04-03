**Tree Convert**

Tree Convert is a tool for converting ISOGG to YFull nomenclature, and vice versa.
For separate comparisons you can use the Streamlit GUI (main.py) which also runs at https://treeconvert.streamlit.app/.
For batch operations you can use treecompare.py. This only requires one input file, a CSV file with three columns: Sample, ISOGG_hg, YFull_hg.

**Installing**

Clone this repository and create a virtual environment. Install the packages in requirements.txt.

**Run Streamlit locally**

If you want to run the Streamlit GUI locally, just run the following command inside your virtual environment:
`streamlit run main.py`

**Run batch tool**

To run the batch function, just run the following command inside your virtual environment:
`python treecompare.py --input input_file.csv --output output_file.csv`
Replace `input_file.csv` and `output_file.csv` with your files.

The batch tools outputs a csv file with the following columns:
- sample: Sample name
- isogg: Input ISOGG haplogroup
- yfull: Input YFull haplogroup
- isogg_path: Full path for ISOGG haplogroup
- yfull_path: Full path for YFull haplogroup
- isogg_identifying_snps: Identifying SNPs for ISOGG haplogroup
- yfull_identifying_snps: Identifying SNPs for YFull haplogroup
- isogg_missing_snps: Identifying SNPs for ISOGG haplogroup that are missing in YFull tree
- yfull_missing_snps: Identifying SNPs for YFull haplogroup that are missing in ISOGG tree
- isogg_hg_yfull_snps: Inferred ISOGG haplogroups based on YFulls identifying SNPs
- yfull_hg_isogg_snps: Inferred YFull haplogroups based on ISOGGs identifying SNPs
- highest_isogg_hg: Consensus ISOGG haplogroup that is supported by the most identifying SNPs
- highest_isogg_hg_ratio: If different identifying SNPs lead to different haplogroups, this ratio is the percentage of SNPs that lead to this consensus highest haplogroup
- highest_yfull_hg: Consensus YFull haplogroup that is supported by the most identifying SNPs
- highest_yfull_hg_ratio: If different identifying SNPs lead to different haplogroups, this ratio is the percentage of SNPs that lead to this consensus highest haplogroup
- isogg_resolution: Compares the consensus ISOGG haplogroup - which is based on the YFull identifying SNPs - to the true ISOGG haplogroup
- yfull_resolution: Compares the consensus YFull haplogroup - which is based on the ISOGG identifying SNPs - to the true YFull haplogroup
- log: List of warnings

**How does it work?**

The conversion is basically based on six JSON-files: `isogg_tree`, `yfull_tree`, `isogg_hg_to_snp_dict`, `yfull_hg_to_snp_dict`, `isogg_snp_to_hg_dict`, and `yfull_snp_to_hg_dict`.
When you enter an ISOGG/Yfull haplogroup it uses the `isogg_tree` and `yfull_tree` to reconstruct the path to the root (Y-Chromosome 'Adam').
Next it uses `isogg_hg_to_snp_dict` and `yfull_hg_to_snp_dict` to get the identifying SNPs for these haplogroups.
Lastly it uses `isogg_snp_to_hg_dict` and `yfull_snp_to_hg_dict` to convert the identifying SNPs back to a haplogroup in _the other tree_.
It then indicated if one of the two trees has a higher resolution, if they have the same resolution, or if the paths do not correspond.

**FAQ**
__Why do isogg_resolution and yfull_resolution not give the same results?__
Often, the isogg_resolution and yfull_resolution do not give the same result. This is expected. For instance, the yfull_resolution says the resolutions are the same, while the isogg_resolution indicates that YFull resolution is 1 level higher. This only means that the inferred YFull haplogroup - based on the ISOGG identifying SNPs - is the same as the true YFull haplogroup, but that the inferred ISOGG haplogroup - based on the YFull identifying SNPs - has one level higher resolution than the true ISOGG haplogroup. It is therefore important to check both columns for the true resolution difference.

Please note that haplogroups ending in ~ (tilde) are considered one level higher resolution than the same haplogroup without ~.
Also always check the log column to see if there are ambiguous consensus haplogroups. If so, it is best to check these comparisons manually using the Streamlit GUI, because this gives you more information.
