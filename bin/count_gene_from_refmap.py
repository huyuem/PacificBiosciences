'''
This script count the number of reads match each gene from the refmap file
given by cuffcompare/gffcompare.

example line:
ref_gene_id	ref_id	class_code	qry_id_list
SPSB1	TCONS_00000144	c	m160614_013904_42139_c100966090630000001823217407061600_s1_p0/95347/1231_65_CCS|m160614_013904_42139_c100966090630000001823217407061600_s1_p0/95347/1231_65_CCS.mrna1
KAZN	TCONS_00000240	c	m160614_141925_42139_c100966090630000001823217407061602_s1_p0/104835/1067_69_CCS|m160614_141925_42139_c100966090630000001823217407061602_s1_p0/104835/1067_69_CCS.mrna1
KAZN	TCONS_00000242	c	m160614_141925_42139_c100966090630000001823217407061602_s1_p0/104835/1067_69_CCS|m160614_141925_42139_c100966090630000001823217407061602_s1_p0/104835/1067_69_CCS.mrna1
KAZN	TCONS_00000243	c	m160614_141925_42139_c100966090630000001823217407061602_s1_p0/104835/1067_69_CCS|m160614_141925_42139_c100966090630000001823217407061602_s1_p0/104835/1067_69_CCS.mrna1

For read matching n genes (n > 1), each of those gene gets 1/n count
only specific type of cuff/gffcompare is allowed, by default '=' and 'c'. It
does exactly the same thing as the following awk command:

gawk 'BEGIN{FS=OFS="\t"}{if(ARGIND==1) {n=split($4,a,","); for(i=1;i<=n;++i) loci[a[i]]++;} else {n=split($4,a,","); for(i=1;i<=n;++i){iso[$2]+=1.0/loci[a[i]]}}}END{for(i in iso) printf "%s\t%.2f\n", i, iso[i]; for(r in loci) printf "%s\t%d\n", r, loci[r] > "/dev/stderr";  }' <(tail -n +2 refmap) <(tail -n +2 refmap) 1> counts 2>loci

The awk on some Linux system does not support variable ARGIND, thus this python
script is provided.

author: Bo Han (bhan@pacb.com)
last modified: 2016-06-30 14:48
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
        for read in tokens[3].split(','):
            counter[read] += 1
            filedata.append([tokens[0].strip(), tokens[1], tokens[2], read])
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
