from rds2py import read_rds

# Load a pre-processed R object (e.g., from a colleague)
r_obj = read_rds("processed_data.rds")

# Uses BiocPy's SingleCellExperiment representation
print(r_obj)

