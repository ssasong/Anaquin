#
# python3 scripts/hg382chrQ.py data/genome/hg38/sequin_regions_hg38_2.6.bed data/genome/chrQ/sequin_regions_chrQ_2.7.bed <Bed>
#

import sys

if len(sys.argv) != 4:
    print("\nUsage: python hg382chrQ.py <Sequin_HG38_BED> <Sequin_ChrQ_BED> <Input BED in HG38>\n\nThe script will convert between hg38 to chrQ coordinates. It is important to note the hg38 coordinates must overlap with <Sequin_HG38_BED> BED file.\nCopy and paste the commands below for a sample demo. Note the input BED file is in hg38 but the output will be in chrQ.\n\nwget https://raw.githubusercontent.com/sequinstandards/Anaquin/master/data/genome/hg38/sequin_regions_hg38_2.6.bed\nwget https://raw.githubusercontent.com/sequinstandards/Anaquin/master/data/genome/chrQ/sequin_regions_chrQ_2.7.bed\nwget https://raw.githubusercontent.com/sequinstandards/Anaquin/master/test/hg38.bed\npython hg382chrQ.py sequin_regions_hg38_2.6.bed sequin_regions_chrQ_2.7.bed hg38.bed\n")
else:
    hg38 = sys.argv[1]
    chrQ = sys.argv[2]
    file = sys.argv[3]

    with open(hg38) as r:
        hg38 = {}
        for line in r:
            toks = line.strip().split("\t")
            hg38[toks[3]] = { "chr":toks[0], "start":int(toks[1]), "end":int(toks[2]), "name":toks[3] }
        
    with open(chrQ) as r:
        chrQ = {}
        for line in r:
            toks = line.strip().split("\t")
            assert(not toks[3] in chrQ)
            chrQ[toks[3]] = { "chr":toks[0], "start":int(toks[1]), "end":int(toks[2]), "name":toks[3] }

    assert(len(hg38) > 0)
    assert(len(chrQ) > 0)
    n = 0
    last = None

    with open(file) as r:
        keys = list(reversed(sorted(hg38.keys())))
        for line in r:
            toks = line.strip().split("\t")
        
            x0 = toks[0]
            x1 = int(toks[1])
            x2 = int(toks[2])
            found = False
            
            if "chrQ" in x0:
                print(line)
                continue

            # It's slow but it's alright...
            for i in keys:
                # Only if it's intersecting, otherwise not possible to convert
                if hg38[i]["chr"] == x0 and hg38[i]["start"] <= x1 and hg38[i]["end"] >= x2:
                    name = hg38[i]["name"].replace("_R", "").replace("_A", "") # Sequin name
                    if not name in chrQ:
                        raise Exception(name + " is in hg38 BED file, but not in chrQ BED file")

                    d1_38 = x1 - hg38[i]["start"] # Delta of the starting position on hg38
                    d2_38 = hg38[i]["end"] - x2   # Delta of the ending position on hg38
                
                    assert(d1_38 >= 0)
                    assert(d2_38 >= 0)

                    d1 = chrQ[name]["start"] + d2_38
                    d2 = chrQ[name]["end"] - d1_38

                    if d1 > d2:
                        raise Exception("Conversion failed. Please contact us for further details")
                    
                    if last != toks[3]:
                        n = 0 # Reset it back
                    last = toks[3]
                    toks[0] = chrQ[name]["chr"]
                    toks[1] = str(d1)
                    toks[2] = str(d2)
                    toks[3] = toks[3] + "_" + str(n)
                    n += 1
                    
                    print('\t'.join(toks))
                    found = True
                    break
        
            if not found:
                raise Exception("Failed to find intersection for: [" + line.replace("\t", " ").strip() + "]")
#<<@@@@>>