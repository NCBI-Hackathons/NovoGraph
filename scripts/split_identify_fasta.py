#!/usr/bin/env python3

from argparse import ArgumentParser
import sys
import re

parser = ArgumentParser()
parser.add_argument("inputfile", help = "Input File in fasta format")
parser.add_argument("prefix", help = "This will be prefixed to the fasta ids")
parser.add_argument("outputfile", help="Output with split fasta records")
parser.add_argument("--Ns", type = int, default = 5000, help = "Number of Ns to split by [5000]")

args = parser.parse_args()

def read_fasta(fp):
    """Parses Fasta file using a generator function."""
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))
     
def get_new_id(prefix,partnr,oldid):
    split_oldid = oldid.split(" ")
    if len(split_oldid) > 1:
        fid = split_oldid[0]
        rid = " ".join(split_oldid[1:])
    else:
        fid = split_oldid[0]
        rid = ""
    return ">" + args.prefix + "_" + fid.lstrip(">") + "_" + str(partnr) + "_" + rid + "\n"

with open(args.outputfile, "w+") as outf:
    with open(args.inputfile) as f:
        for fid, seqt in read_fasta(f):
            seq = seqt.upper()
            if seq[0] not in "GATCRYWSMKHBVDN":
                print("Error found in input file. Make sure the input is in correct fasta format.")
                sys.exit()
            print(fid + "... ", end = "")

            s = seq.rstrip().upper().lstrip("N").rstrip("N")
            matches = re.finditer("N{" + str(args.Ns) + ",}", s)
            start = 0
            partnr = 1
            for m in matches:
                nstart, nstop = m.span()
                #print(s[m.span()[1]])
                part = s[start:nstart]

                outf.write(get_new_id(args.prefix, str(partnr),fid))
                outf.write(part + "\n")

                partnr += 1
                start = nstop
            part = s[start:]
            outf.write(get_new_id(args.prefix, str(partnr),fid))
            outf.write(part + "\n")
            
            print("found " + str(partnr) + " parts")

