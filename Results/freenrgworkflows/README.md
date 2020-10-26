### Input/output to freenrgworkflows network analysis

See script for usage, eg:

`python run_networkanalysis.py group1/G1_AMBER_summary.csv --target_compound Lig5 -o out.dat -e group1/Group1_ic50_exp.csv --stats --generate_notebook`

takes the raw experimental and computed AMBER data for group 1, and computes DDGs and their errors, all relative to Lig5.
