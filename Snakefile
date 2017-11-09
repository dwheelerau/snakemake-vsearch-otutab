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
CLUSTERMODE = config['CLUSTERMODE']
CLUSTERID = config['CLUSTERID']
REF_DB = config['REF_DB']
TAX_REF = config['TAX_REF']

rule all:
    input: "sample.otu_table.txt", "tax_summary_proportion/"


#rule project_setup:
#    output: DIRS
#    shell:
#        "mkdir -p "+' '.join(DIRS)

rule cluster_otus:
    input:
        "scripts/run_otu.sh"
    output:
        "sample.otu_table.txt"
    log:
        stdout=expand("run_otu.{date_string}.stdout",date_string=date_string),
        stderr=expand("run_otu.{date_string}.stderr",date_string=date_string)
    params:
        threads=expand("{THREADS}", THREADS=THREADS),
        minovlen=expand("{MINOVLEN}", MINOVLEN=MINOVLEN),
        maxdiff=expand("{MAXDIFF}", MAXDIFF=MAXDIFF),
        maxee=expand("{MAXEE}", MAXEE=MAXEE),
        minlen=expand("{MINLEN}", MINLEN=MINLEN),
        maxlen=expand("{MAXLEN}", MAXLEN=MAXLEN),
        maxns=expand("{MAXNS}", MAXNS=MAXNS),
        clustermode=expand("{CLUSTERMODE}", CLUSTERMODE=CLUSTERMODE),
        clusterid=expand("{CLUSTERID}", CLUSTERID=CLUSTERID),
        ref_db=expand("{REF_DB}", REF_DB=REF_DB)
    shell:
        "scripts/run_otu.sh "
        "-t {params.threads} -m {params.minovlen} -d {params.maxdiff} "
        "-e {params.maxee} -i {params.minlen} -l {params.maxlen} "
        "-n {params.maxns} -c '{params.clustermode}' -r {params.clusterid} "
        "-b '{params.ref_db}' > >(tee -a {log.stdout}) 2> >(tee -a {log.stderr} >&2)"

rule conv_biom:
    input:
        "sample.otu_table.txt"
    output:
        "sample.otu_table.biom"
    shell:
        "biom convert -i {input} -o {output} --table-type 'OTU table' --to-hdf5"

rule assign_taxonomy:
    input:
        "all.otus.fasta" 
    output:
        "taxonomy/all.otus_tax_assignments.txt"
    params:
        taxonomy=expand("{TAX_REF}", TAX_REF=TAX_REF),
        taxref=expand("{REF_DB}", REF_DB=REF_DB)
    shell:
        "assign_taxonomy.py -m uclust -i {input} -o taxonomy "
        "-t {params.taxonomy} -r {params.taxref}"

rule add_taxonomy:
    input:
        taxa="taxonomy/all.otus_tax_assignments.txt",
        biom="sample.otu_table.biom"
    output:
        "sample.otu_table_wTax.biom"
    shell:
        "biom add-metadata -i {input.biom} --observation-metadata-fp {input.taxa} "
        "--sc-separated taxonomy --observation-header OTUID,taxonomy -o {output}"

rule summarize_taxa:
    input:
        "sample.otu_table_wTax.biom"
    output:
        prop="tax_summary_proportion/",
        count="tax_summary_count/"
    shell:
        """
        summarize_taxa.py -i {input} -o {output.prop} --suppress_biom_table_output
        summarize_taxa.py -a -i {input} -o {output.count} --suppress_biom_table_output
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
        """
