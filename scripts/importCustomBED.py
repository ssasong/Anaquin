import os
import sys

with open(sys.argv[1]) as r:
    for line in r:
        toks = line.strip().split("\t")
        toks[3] = toks[3] + "_" + toks[1] + "_" + toks[2]
        print("\t".join(toks))
