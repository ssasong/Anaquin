#
# python3 scripts/unflip.py <Input> <Output>
#

import sys
import pysam
from Bio import SeqIO

def comp(s): 
    basecomplement = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N'} 
    letters = list(s) 
    letters = [basecomplement[base] for base in letters] 
    return ''.join(letters)

def rComp(s):
    return comp(s[::-1])

def reverse(s):
    return s[::-1]

if sys.argv[1].endswith(".bam"):
    r = pysam.AlignmentFile(sys.argv[1], "rc")
    w = pysam.AlignmentFile(sys.argv[2], "wb", template=r)
    for i in r:
        #i.seq = i.seq[::-1]
        i.seq = comp(i.seq)
        if i.is_reverse:
            i.seq = rComp(i.seq)
        w.write(i)
    r.close()
    w.close()
else:
    with open(sys.argv[2], "w") as w:
        with open(sys.argv[1], "r") as r:
            i = 0
            for line in r:
                if i == 1:
                    w.write(comp(line.strip()) + "\n")
                else:
                    w.write(line.strip() + "\n")
                i = (i + 1) % 4
