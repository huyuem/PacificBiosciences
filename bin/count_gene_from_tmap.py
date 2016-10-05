'''
This script count the number of reads match each gene from the tmap file
given by cuffcompare/gffcompare.

example line:
TMCO1	NM_019026	=	m160403_065056_42175_c100993391270000001823222007191686_s1_p0/50/1531_55_CCS	m160403_065056_42175_c100993391270000001823222007191686_s1_p0/50/1531_55_CCS.mrna1	10	0.000000	0.000000	0.000000	0.000000	1491	m160403_065056_42175_c100993391270000001823222007191686_s1_p0/50/1531_55_CCS.mrna1	4470

For read matching n genes (n > 1), each of those gene gets 1/n count
only specific type of cuff/gffcompare is allowed, by default '=' and 'c'. It
does exactly the same thing as the following awk command:

gawk 'BEGIN{FS=OFS="\t"}{if(ARGIND==1)a[$4]++; else{ if($3=="=" || $3=="c") c[$2]+=1/a[$4]}}END{for(g in c) print g, c[g]}' file file

The awk on some Linux system does not support variable ARGIND, thus this python
script is provided.

author: Bo Han (bhan@pacb.com)
last modified: 2016-04-10 10:10
'''
import sys
from collections import Counter, defaultdict

def count_read_hits(file):
    counter = Counter()
    fh = open(file)
    fh.readline() # consume the header
    filedata = []
    for line in fh.readlines():
        tokens = line.split('\t')
        counter[tokens[3]] += 1
        filedata.append([tokens[0].strip(), tokens[1], tokens[2], tokens[3] ])
    return counter, filedata # can afford this in memory

def count_mRNA_hits(filedata, reads_count, allowed_type = ("=", "c")):
    gene_counter = defaultdict(float)
    for tokens in filedata:
        if tokens[2] in allowed_type:
            gene_counter[tokens[1]] += 1.0/reads_count[tokens[3]]
    return gene_counter

def count_gene_hits(filedata, reads_count, allowed_type = ("i", "j", "x", "o", "p", "c", "=", "e")):
    gene_counter = defaultdict(float)
    for tokens in filedata:
        if tokens[2] in allowed_type:\
            gene_counter[tokens[0]] += 1.0/reads_count[tokens[3]]
    return gene_counter

def report(ct):
    for g in ct:
        print("{}\t{}".format(g, ct[g]))

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("usage: {} cuff/gffcompare.tmap gene|mRNA".format(sys.argv[0]))
    if sys.argv[2] == "mRNA":
        reads_counter, fd = count_read_hits(sys.argv[1])
        genes_counter = count_mRNA_hits(fd, reads_counter)
    elif sys.argv[2] == "gene":
        reads_counter, fd = count_read_hits(sys.argv[1])
        genes_counter = count_gene_hits(fd, reads_counter)
    else:
        print("unsupported class, use only gene or mRNA")
        sys.exit(1)
    report(genes_counter)
