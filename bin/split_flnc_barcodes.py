from Bio import SeqIO
import sys

if len(sys.argv) < 2:
    sys.stderr.write("usage: {} flnc.fa num_barcodes[8]".format(sys.argv[0]))
    sys.exit(1)

if len(sys.argv) == 3:
    num_bar = int(sys.argv[2])
else:
    num_bar = 8
    
fhs = {}
for n in xrange(0, num_bar):
    fhs[str(n)] = open("{}.barcode{}".format(sys.argv[1], n), 'w')

for fa in SeqIO.parse(open(sys.argv[1]), "fasta"):
    primerinfo = fa.description.split(' ')[1].split(";")[7].split('=')
    assert primerinfo[0] == 'primer'
    if primerinfo[1] != "NA":
        SeqIO.write(fa, fhs[primerinfo[1]], "fasta")