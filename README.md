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

**How does it work?**

The conversion is basically based on six JSON-files: `isogg_tree`, `yfull_tree`, `isogg_hg_to_snp_dict`, `yfull_hg_to_snp_dict`, `isogg_snp_to_hg_dict`, and `yfull_snp_to_hg_dict`.
When you enter an ISOGG/Yfull haplogroup it uses the `isogg_tree` and `yfull_tree` to reconstruct the path to the root (Y-Chromosome 'Adam').
Next it uses `isogg_hg_to_snp_dict` and `yfull_hg_to_snp_dict` to get the identifying SNPs for these haplogroups.
Lastly it uses `isogg_snp_to_hg_dict` and `yfull_snp_to_hg_dict` to convert the identifying SNPs back to a haplogroup in _the other tree_.
