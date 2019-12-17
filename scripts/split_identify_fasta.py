from argparse import ArgumentParser
import sys
import re

parser = ArgumentParser()
parser.add_argument("inputfile", help = "Input File in fasta format")
parser.add_argument("prefix", help = "This will be prefixed to the fasta ids")
parser.add_argument("outputfile", help="Output with split fasta records")
parser.add_argument("--Ns", type = int, default = 5000, help = "Number of Ns to split by [5000]")

args = parser.parse_args()

with open(args.outputfile, "w+") as outf:
    with open(args.inputfile) as f:
        for linenr, line in enumerate(f):
            if linenr % 2 == 0:
                rid = line.rstrip().lstrip(">")
            else:
                if line[0] not in {'A', 'G', 'C', 'T', 'N'}:
                    print("Error found in input file. Make sure the input is in correct fasta format.")
                    sys.exit()
                print(rid + "... ", end = "")

                s = line.rstrip().upper().lstrip("N").rstrip("N")
                matches = re.finditer("N{" + str(args.Ns) + ",}", s)
                start = 0
                partnr = 1
                for m in matches:
                    nstart, nstop = m.span()
                    #print(s[m.span()[1]])
                    part = s[start:nstart]

                    outf.write(">" + args.prefix + "_" + rid + "_" + str(partnr) + "\n")
                    outf.write(part + "\n")

                    partnr += 1
                    start = nstop
                part = s[start:]
                outf.write(">" + args.prefix + "_" + rid + "_" + str(partnr) + "\n")
                outf.write(part + "\n")
                
                print("found " + str(partnr) + " parts")

