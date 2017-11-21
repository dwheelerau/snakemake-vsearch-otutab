#!/usr/bin/env python2
import sys

# simple script to rename merged sample names into something better
# requires a sample_key.tsv file that looks like this (short name \t long name)
# A40C    A-40-C
# A41D    A-41-D
# the script will use periods which are allowed by qiime
# ie A40C -> A.40.C

try:
    fname = sys.argv[1]
    key_dict = dict(((line.split('\t')[0], line.strip().split('\t')[1])
                     for line in open(fname)))
except:
    print('Missing keyfile?')
    sys.exit(1)

outfile = open('sample.rename.otu_table.txt', 'w')

header = []
with open('sample.otu_table.txt') as f:
    for line in f:
        if len(header) == 0:
            new_header = []
            header = line.strip().split('\t')
            for name in header:
                if len(new_header) == 0:
                    new_header.append(name)  # OTU ID
                else:
                    try:
                        new_name = key_dict[name]
                        new_name = new_name.replace('-', '.')
                        new_header.append(new_name)
                    except KeyError:
                        print('key file doesnt contain all samples name pairs')
                        print('current key file is: %s' % fname)
                        print line
                        sys.exit(1)
            new_header = '\t'.join(new_header) + '\n'
            outfile.write(new_header)
        else:
            outfile.write(line)
outfile.close()
