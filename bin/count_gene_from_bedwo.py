'''
This script count the number of reads match each gene from the bedwo file
given by bedtools intersect -wo.

example line:
gi|48994873|gb|U00096.2|	336	2799	b0002	0	+	336	2799	0	1	2463,	0,	gi|48994873|gb|U00096.2|	301	1881	m160304_022746_42175_c100978430030000001823222908031681_s1_p0/4515/1641_60_CCS	60	+	1545

For read matching n genes (n > 1), each of those gene gets 1/n count

author: Bo Han (bhan@pacb.com)
last modified: 2016-04-15 15:22
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
        counter[tokens[15]] += 1
        filedata.append([tokens[3], tokens[5], tokens[15], tokens[17] ])
    return counter, filedata # can afford this in memory

class count(object):
    __slots__ = ("sense_count", "antisense_count")

    def __init__(self, sc = 0.0, asc = 0.0):
        self.sense_count = sc
        self.antisense_count = asc

    def add_sense(self, c):
        self.sense_count +=c

    def add_antisense(self, c):
        self.antisense_count += c

def count_gene_hits(filedata, reads_count):
    ct = defaultdict(count)
    for tokens in filedata:
        if tokens[1] == tokens[3]:
            ct[tokens[0]].add_sense(1.0/reads_count[tokens[2]])
        else:
            ct[tokens[0]].add_antisense(1.0/reads_count[tokens[2]])
    return ct

def report(ct):
    for g in ct:
        print("{}\t{}\t{}".format(g, ct[g].sense_count, ct[g].antisense_count))

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("usage: {} bedwo".format(sys.argv[0]))
    reads_counter, fd = count_read_hits(sys.argv[1])
    genes_counter = count_gene_hits(fd, reads_counter)
    report(genes_counter)
