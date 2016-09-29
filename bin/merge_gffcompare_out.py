'''
this script organizes the output of gffcompare and generate a table for R to plot
The input looks like:
------------------------------------------------------------------------------------
# gffcompare v0.9.6 | Command line was:
# gffcompare  ...
#

#= Summary for dataset: ../sampleA.all_sizes.quivered_hq.fastq.5merge.collapsed.min_fl_2.gff
#     Query mRNAs :    7995 in    6247 loci  (6447 multi-exon transcripts)
#            (1092 multi-transcript loci, ~1.3 transcripts per locus)
# Reference mRNAs :   12027 in    5607 loci  (11563 multi-exon)
# Super-loci w/ reference transcripts:     5257
#-----------------| Sensitivity | Precision  |
        Base level:    70.3     |    66.5    |
        Exon level:    69.0     |    87.2    |
      Intron level:    70.0     |    96.3    |
Intron chain level:    33.8     |    60.5    |
  Transcript level:    32.9     |    49.5    |
       Locus level:    65.8     |    59.3    |

     Matching intron chains:    3903
       Matching transcripts:    3956
              Matching loci:    3689

          Missed exons:   16534/72074	( 22.9%)
           Novel exons:    1928/58000	(  3.3%)
        Missed introns:   14345/68091	( 21.1%)
         Novel introns:     372/49481	(  0.8%)
           Missed loci:       0/5607	(  0.0%)
            Novel loci:     370/6247	(  5.9%)

#= Summary for dataset: ../sampleB.all_sizes.quivered_hq.fastq.5merge.collapsed.min_fl_2.gff
#     Query mRNAs :    6545 in    5151 loci  (5641 multi-exon transcripts)
#            (913 multi-transcript loci, ~1.3 transcripts per locus)
# Reference mRNAs :   10787 in    4975 loci  (10387 multi-exon)
# Super-loci w/ reference transcripts:     4701
#-----------------| Sensitivity | Precision  |
        Base level:    69.3     |    79.6    |
        Exon level:    69.5     |    88.5    |
      Intron level:    70.3     |    96.7    |
Intron chain level:    35.4     |    65.1    |
  Transcript level:    34.6     |    57.1    |
       Locus level:    69.9     |    67.7    |

     Matching intron chains:    3673
       Matching transcripts:    3736
              Matching loci:    3476

          Missed exons:   13814/61908	( 22.3%)
           Novel exons:    1038/49431	(  2.1%)
        Missed introns:   11815/58429	( 20.2%)
         Novel introns:     226/42473	(  0.5%)
           Missed loci:       0/4975	(  0.0%)
            Novel loci:     138/5151	(  2.7%)

 Total union super-loci across all input datasets: 7628
  (1796 multi-transcript, ~1.5 transcripts per locus)
'''

import sys
from os.path import basename

class StatObj(object):
    '''wrap sensitivity and precision'''
    __all__ = ("sensitivity", "precision")
    def __init__(self, (s, p,)):
        try:
            self.sensitivity = float(s)
            self.precision = float(p)
        except:
            print("cannot parse {} and {}".format(s, p))

    def __str__(self):
        '''for easy print'''
        return str(self.sensitivity) + '\t' + str(self.precision)
    __repr__ = __str__

    def __getitem__(self, attr):
        return self.__dict__[attr]

class StatSum(object):
    '''parsed statistical data for each sample'''
    __all__ = ("sample_name", "base", "exon", "intron", "intron_chain", "transcript", "locus")

    def __init__(self, fh):
        l = fh.readline()
        if not l or l[0] != '#':
            raise EOFError('')
        self.sample_name = basename(l.strip().split()[4])
        for i in range(0, 5): # skip the next 5 lines for this entry
            fh.readline()
        self.base = StatObj(self._from_fh(fh))
        self.exon = StatObj(self._from_fh(fh))
        self.intron = StatObj(self._from_fh(fh))
        self.intron_chain = StatObj(self._from_fh(fh, 3, 5)) # intron chain is split into two tokens
        self.transcript = StatObj(self._from_fh(fh))
        self.locus = StatObj(self._from_fh(fh))
        for i in range(0, 12): # skip the next 6 lines
            fh.readline()

    @staticmethod
    def _from_fh(fh, i = 2, j = 4):
        tokens = fh.readline().strip().split()
        return tokens[i], tokens[j]

    def __getitem__(self, attr):
        return self.__dict__[attr]

def parse(file):
    outs = []
    fh = open(file)
    for i in range(0, 4):
        fh.readline() # skip the header of the file
    while True:
        try:
            s = StatSum(fh)
        except EOFError:
            break
        outs.append(s)
    return outs

def report(lst):
    print("category\tparameter\t"),
    print('\t'.join([s.sample_name for s in lst]))
    for attr in StatSum.__all__[1:]:
        for attr2 in StatObj.__all__:
            print("{}\t{}\t".format(attr, attr2)),
            print('\t'.join([str(s[attr][attr2]) for s in lst]))

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("usage:\n\t{} combined.gffcompare".format(sys.argv[0]))
        sys.exit(1)
    report(parse(sys.argv[1]))
