import scanpy as sc

adata = sc.read_10x_mtx(
    '../inferCSN_data/science_data/rna/',  
    var_names='gene_symbols',                
    cache=True)
