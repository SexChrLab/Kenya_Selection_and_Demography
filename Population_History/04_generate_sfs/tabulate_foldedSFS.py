import argparse
import os

parser = argparse.ArgumentParser(description="Combine results of sfs from all of the chromosomes on the autosomes.")
parser.add_argument("--num_bin",required=True,help="Input the number of bins of the sfs.")
parser.add_argument("--directory",required=True,help="Input the directory where the sfs files are.")
parser.add_argument("--out_filename",required=True,help="Input the name of the output file.")

args = parser.parse_args()

num=int(args.num_bin)
autosome_eta = {}

for i in range(1, num+1):
    autosome_eta[i] = 0

for i in range(1,23):
    filename = os.path.join(args.directory, "chr" + str(i) + "_sfs.txt")
    with open(filename, "r") as f:
        for line in f:
            if not line.startswith("af_bin"):
                line = line.rstrip("\n").split("\t")
                autosome_eta[int(line[0])] += float(line[1])

outfile = open(args.out_filename, "w")
for k, v in autosome_eta.items():
    toprint = [str(k), str(v)]
    print ("\t".join(toprint), file=outfile)
