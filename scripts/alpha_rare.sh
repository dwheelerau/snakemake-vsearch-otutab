# -e sets level, random sampling to this level across all samples, any samples
# with less than this value are dropped from the table.
alpha_rarefaction.py -i ./otu_table_newIDs_w_tax.filtered.biom -o arare_max5000/ -t ./rep_set_numbered.tre -m NZGL02039_meta_fixed2.tsv -e 5000 
