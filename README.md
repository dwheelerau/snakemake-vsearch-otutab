## Snakemake pipeline for vsearch based qiime analysis  


### Requirments  
Paired end reads located in a directory called `projectData` located in the same
path as the Snakefile from this repo. The filenames need to be in the following
form:

`SAMPLENAME_L001_R1_001.fastq.gz`

The SAMPLENAME can not contain ANY special characters as the script uses
periods (.) and underscores to parse out the sample name information away from
the rest of the file name.  
-  snakemake  
-  vsearch  
-  QIIME1 (I have this installed via conda)  

The final steps of the analysis tries to fix the sample names by adding a
period (.) in place of any original mutations that were required to avoid
special characters (periods are allowed via qiime). To do this a python script
called `fix_names.py` in the scripts directory is executed. This file requires
a `sample_key.tsv` file in the following format.  

```
# NOTE: the script will actually use periods instead of dashes  
F101110 F-10-11-10
F101110 F-10-111-10
```

Creation of the biom table also requires a meta file in qiime format, this
should use the period names that will be outputed by the `fix_names.py` script.
The meta file name must be specified in the `config.yaml` file. Here is an example of what this workflow will do to an example sample names.  
  
```
Sample-1 -> sample1 -> sample.1  
# meta file should look like this  

#SampleID    BarcodeSequence     LinkerPrimerSequence    Description  
Sample.1 aaaaaaaaaaaaaaaaaa  aaaaaaaaaaaaaaaaaa  Sample-1  
```  

*Importantly* you need to adjust the config file to meet your requirements.  

### Running the snakemake pipeline
`snakemake all` will run the full pipeline.  
`snakemake clean` will remove all outputs.

Specify the number of cps you want to expose to the pipeline using `--cores`.  
`snakemake all --cores 24`  


### The rules available  
Individual rule can be run with `snakemake *rule_name*`  

`rule all` - run all rules.  
`rule cluster_otus` - cluster otus using vsearch.  
`rule conv_biom` - convert OTU table to biom table.  
`rule assign_taxonomy` - assign taxonomy using green genes.  
`rule add_taxonomy` - add taxonomy to biom table.  
`rule summarize_taxa` - summarise taxa information.  
`rule summarize_table` 
`rule align_otus` - align OTU sequences for phylogeny.  
`rule filter_aln`  
`rule make_phylogeny` - make phylogeny from alignment.  
`rule core_div_analysis` - run core_diversity_analysis.py script.  
`rule clean` - remove all files for fresh run.  


### testing  
Information on Mock communities can be downloaded from here:
`https://raw.githubusercontent.com/caporaso-lab/mockrobiota/master/inventory.tsv`  


### Config.yaml file  
The config file is shown below, note important options for overlap and also key
params for the `core_diversity_analysis.py` script in qiime. The snakemake
pipeline will output the sample counts for determining the `SAM_DEPTH` value.

```
## general params
REF_DB: /home/dwheeler/databases/gg_13_8_otus/rep_set/97_otus.fasta
TAX_REF: /home/dwheeler/databases/gg_13_8_otus/taxonomy/97_otu_taxonomy.txt
THREADS: 24

## MERGE PARAMs
# minimum overlap between pairs
MINOVLEN: 60
# max differences in overlap
MAXDIFF: 3

## FILTER PARAMs
# max errors
MAXEE: 0.5
# min length of resulting merged fragment
MINLEN: 400
# max length of resulting merged fragment
MAXLEN: 500
# max number of Ns
MAXNS: 0

## CLUSTER PARAMs
CLUSTERMODE: --cluster_size
CLUSTERID: 0.97

## Core diversity analysis settings, using a default of 5000
SAM_DEPTH: 100
MAP_FILE: meta.tsv
```
