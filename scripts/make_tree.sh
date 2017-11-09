# activate
#qiime1

# Align sequences command 
#align_seqs.py -i ../rep_set_numbered.fa -o pynast_aligned_seqs/

# Filter alignment command 
#filter_alignment.py -o pynast_aligned_seqs/ -i ./pynast_aligned_seqs/rep_set_numbered_aligned.fasta

# Build phylogenetic tree command 
make_phylogeny.py -i ./pynast_aligned_seqs/rep_set_numbered_aligned_pfiltered.fasta \
  -o rep_set_numbered.tre
