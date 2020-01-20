
from argparse import ArgumentParser
import sys
import re

parser = ArgumentParser()
parser.add_argument("inputfile", help = "Input File in vcf format")
parser.add_argument("contig_alignment_file", help = "Alignments per chromosome")
parser.add_argument("id_folder", help = "Directory containing *.ids")
parser.add_argument("outputfile", help="Output 

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

                outf.write(">" + args.prefix + "_" + fid.lstrip(">") + "_" + str(partnr) + "\n")
                outf.write(part + "\n")

                partnr += 1
                start = nstop
            part = s[start:]
            outf.write(">" + args.prefix + "_" + fid.lstrip(">") + "_" + str(partnr) + "\n")
            outf.write(part + "\n")
            
            print("found " + str(partnr) + " parts")

