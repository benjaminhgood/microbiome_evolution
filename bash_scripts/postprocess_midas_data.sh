species=$1

python ~/projectBenNandita/calculate_marker_gene_coverage.py $species
#python ~/projectBenNandita/combine_replicates.py $species
python ~/projectBenNandita/calculate_coverage_distribution.py $species
python ~/projectBenNandita/calculate_error_pvalues.py $species
