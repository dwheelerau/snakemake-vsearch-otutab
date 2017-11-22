from os.path import join
import glob
from datetime import datetime, date, time
configfile: "./config.yaml"

# globals
now = datetime.now().strftime('%Y-%m-%d-%H-%M')
date_string = "{}".format(now)

# from config.yaml file
THREADS = config['THREADS']
MINOVLEN = config['MINOVLEN']
MAXDIFF = config['MAXDIFF']
MAXEE = config['MAXEE']
MINLEN = config['MINLEN']
MAXLEN = config['MAXLEN']
MAXNS = config['MAXNS']
MINQ = config['MINQ'] # new
CLUSTERMODE = config['CLUSTERMODE']
CLUSTERID = config['CLUSTERID']
REF_DB = config['REF_DB']
TAX_REF = config['TAX_REF']
SAM_DEPTH = config['SAM_DEPTH']
MAP_FILE = config['MAP_FILE']
TAX_METH = config['TAX_METH']

rule all:
    input:
        "sample.otu_table.txt",
        "sample.rename.otu_table.txt",
        "fastqc/",
        "tax_summary_proportion/",
        "sample.rename.otu_table_summary.txt",
        "all.otus.tre",
        "coreout/"

rule cluster_otus:
    input:
        "scripts/run_otu.sh"
    output:
        "sample.otu_table.txt", "all.otus.fasta"
    log:
        stdout=expand("run_otu.{date_string}.stdout",date_string=date_string),
        stderr=expand("run_otu.{date_string}.stderr",date_string=date_string)
    params:
        threads=THREADS,
        minovlen=MINOVLEN,
        maxdiff=MAXDIFF,
        maxee=MAXEE,
        minlen=MINLEN,
        maxlen=MAXLEN,
        maxns=MAXNS,
        minq=MINQ,
        clustermode=CLUSTERMODE,
        clusterid=CLUSTERID,
        ref_db=REF_DB
    shell:
        "scripts/run_otu.sh "
        "-t {params.threads} -m {params.minovlen} -d {params.maxdiff} "
        "-e {params.maxee} -i {params.minlen} -l {params.maxlen} "
        "-n {params.maxns} -q {params.minq} -c '{params.clustermode}' -r {params.clusterid} "
        "-b '{params.ref_db}' > >(tee -a {log.stdout}) 2> >(tee -a {log.stderr} >&2)"

rule qc_reads:
    input:
        "sample.otu_table.txt"
    output:
        "fastqc/"
    shell:
        """
        mkdir -p fastqc
        for fq in tmp/*merged.fastq; do
            fastqc -o fastqc/ $fq
        done
        """

rule fix_names:
    input:
        biom="sample.otu_table.txt",
        key="sample_key.tsv"
    output:
        "sample.rename.otu_table.txt"
    shell:
        "python2 scripts/fix_names.py {input.key}"

rule conv_biom:
    input:
        "sample.rename.otu_table.txt"
    output:
        "sample.rename.otu_table.biom"
    shell:
        "biom convert -i {input} -o {output} --table-type 'OTU table' --to-hdf5"

rule assign_taxonomy:
    input:
        "all.otus.fasta"
    output:
        "taxonomy/all.otus_tax_assignments.txt"
    params:
        taxonomy=TAX_REF,
        taxref=REF_DB,
        taxmeth=TAX_METH
    shell:
        "assign_taxonomy.py -m {params.taxmeth} -i {input} -o taxonomy "
        "-t {params.taxonomy} -r {params.taxref}"

rule add_taxonomy:
    input:
        taxa="taxonomy/all.otus_tax_assignments.txt",
        biom="sample.rename.otu_table.biom"
    output:
        "sample.rename.otu_table_wTax.biom"
    shell:
        "biom add-metadata -i {input.biom} --observation-metadata-fp {input.taxa} "
        "--sc-separated taxonomy --observation-header OTUID,taxonomy -o {output}"

rule summarize_taxa:
    input:
        "sample.rename.otu_table_wTax.biom"
    output:
        prop="tax_summary_proportion/",
        count="tax_summary_count/"
    shell:
        """
        summarize_taxa.py -i {input} -o {output.prop} --suppress_biom_table_output
        summarize_taxa.py -a -i {input} -o {output.count} --suppress_biom_table_output
        """

rule summarize_table:
    input:
        "sample.rename.otu_table_wTax.biom"
    output:
        "sample.rename.otu_table_summary.txt"
    shell:
        """
        biom summarize-table -i {input} -o {output}
        cat {output}
        echo
        echo ' ~~~~~~ Caution! ~~~~~~ '
        echo 'you need the above data to set the SAM_DEPTH (--sampling_depth)'
        echo 'param in config for the core_diversity_analysis rule'
        """
# Align sequences command
rule align_otus:
    input:
        "all.otus.fasta"
    output:
        "pynast_aligned_seqs/all.otus_aligned.fasta"
    shell:
        "align_seqs.py -i {input} -o pynast_aligned_seqs/"

rule filter_aln:
    input:
        "pynast_aligned_seqs/all.otus_aligned.fasta"
    output:
        "pynast_aligned_seqs/all.otus_aligned_pfiltered.fasta"
    shell:
        "filter_alignment.py -o pynast_aligned_seqs/ -i {input}"

rule make_phylogeny:
    input:
        "pynast_aligned_seqs/all.otus_aligned_pfiltered.fasta"
    output:
        "all.otus.tre"
    shell:
        "make_phylogeny.py -i {input} -o {output}"

rule core_div_analysis:
    input:
        summary="sample.rename.otu_table_summary.txt",
        biom="sample.rename.otu_table_wTax.biom",
        meta=MAP_FILE,
        tree="all.otus.tre"
    output:
        "coreout/"
    params:
        sample_depth=SAM_DEPTH
    shell:
        """
        core_diversity_analyses.py -i {input.biom} -o tmpcore/ -m {input.meta} \
        -t {input.tree} -e {params.sample_depth}
        mv tmpcore/* {output}
        rmdir tmpcore/
        """

onerror:
        print("Biom error? Have you run qiime1 from the commmand line to activate it?")
        #shell("mail -s "an error occurred" youremail@provider.com < {log}")

rule clean:
    shell:
        """
        rm -f tmp/*
        rm -f logs/*
        rm -f readStats/*
        rm -f *.uc
        rm -f *.fasta
        rm -f *.fastq
        rm -f all*
        rm -f *.txt
        rm -f *.stats
        rm -f *.biom
        rm -f *.log
        rm -rf tax*
        rm -rf pynast_aligned_seqs/
        rm -rf coreout/
        rm -rf fastqc/
        rm -rf tmpcore/
        """
