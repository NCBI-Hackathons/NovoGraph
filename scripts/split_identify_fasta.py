from argparse import ArgumentParser
from Bio import SeqIO
import re

parser = ArgumentParser()
parser.add_argument("inputfile", help = "Input File in fasta format")
parser.add_argument("prefix", help = "This will be prefixed to the fasta ids")
parser.add_argument("outputfile", help="Output with split fasta records")
parser.add_argument("--Ns", type = int, default = 5000, help = "Number of Ns to split by [5000]")

args = parser.parse_args()

with open(args.outputfile, "w+") as outf:
    for read in SeqIO.parse(args.inputfile, "fasta"):
        print(read.id)

        s = str(read.seq.upper().lstrip("N").rstrip("N"))
        matches = re.finditer("N{" + str(args.Ns) + ",}", s)
        start = 0
        partnr = 1
        for m in matches:
            nstart, nstop = m.span()
            #print(s[m.span()[1]])
            part = s[start:nstart]

            outf.write(">" + args.prefix + "_" + read.id + "_" + str(partnr) + "\n")
            outf.write(part + "\n")

            partnr += 1
            start = nstop
        part = s[start:]
        outf.write(">" + args.prefix + "_" + read.id + "_" + str(partnr) + "\n")
        outf.write(part + "\n")
        
        print("found " + str(partnr) + " parts")

