from argparse import ArgumentParser
import sys
from os import path
import re
from collections import defaultdict, namedtuple

parser = ArgumentParser()
parser.add_argument("inputfile", help = "Input file in vcf format. Only works for overlapping vcf.")
parser.add_argument("contig_alignment_file", help = "Alignments per chromosome")
parser.add_argument("id_prefix", help = "Path and prefix to *.queryIDs")
parser.add_argument("pedigree", help = "File listing HaplotypeIDs to GenotypeIds")
parser.add_argument("outputfile", help="Output file in vcf format")

args = parser.parse_args()

Gtpos = namedtuple('Gtpos', ['gid', 'pos'])

class Gtstring():
    def __init__(self, size):
        self.gtarray = ["."] * size # diploid

    def set(self, pos, val):
        self.gtarray[pos] = str(val)


    def __str__(self):
        final = ""
        for i in range(0, len(self.gtarray), 2):
            final += self.gtarray[i]
            final += "|"
            final += self.gtarray[i+1]
            final += "\t"
        return final.rstrip()
    
class Genotypes():
    def __init__(self,gids,pedigreefile):
        self.gids = gids
        self.htpos = {}
        with open(pedigreefile) as f:
            for line in f:
                gid, sid1, sid2 = line.rstrip().split()
                if gid not in self.gids:
                    print(gid + " in pedigree file not in specified ids.")
                    self = None
                pos = self.gids.index(gid)*2
                self.htpos[sid1] = pos
                self.htpos[sid2] = pos + 1
    def header(self):
        return "\t".join(self.gids)
    
    def get_genotypes(self, chrom, info, table):
        #print(info)
        idfields = info.lstrip("CONTIG=").split(",")
        gtstring = Gtstring(len(self.gids)*2)
        for allelenr, ids in enumerate(idfields):
            print("ids: "  + ids)
            if ids  == "-1":
                continue
            for idx in ids.split("/"):
                if ids  == -1:
                    continue
                s = table[chrom][idx].sampleid
                print("htpos: " + str(self.htpos[s]) + " s: " + str(s))
                gtstring.set(self.htpos[s], allelenr)#print(s)
                #print(self.htpos[s])
                #print(table[chrom][idx].pos)
        return str(gtstring)
        

idfiles = {}
with open(args.contig_alignment_file) as f:
    for line in f:
        if line.startswith("referenceID"):
            continue
        chrom = line.split()[0]
        ipath = args.id_prefix + chrom + ".queryIDs"
        idfiles[ipath] = chrom

idPair = namedtuple('idPair', ['fullid', 'sampleid']) # used to save both the full id and the sampleID

with open(args.pedigree) as f:
    gts = []
    for line in f:
        gt = line.rstrip().split()[0]
        gts.append(gt)

gto = Genotypes(gts, args.pedigree)

#print(gto.dgtpos)

ids_per_chrom = defaultdict(dict)
for idf, chrom in idfiles.items():
    if not path.exists(idf):
        print("Could not find " + str(idf))
        continue
    else:
        with open(idf) as f:
            for line in f:
                fid, vid = line.rstrip().split()
                ids_per_chrom[chrom][vid] = idPair(fid, fid.split("_")[0])


header = ""
with open(args.outputfile, 'w+') as outf:
    with open(args.inputfile) as f:
        for line in f:
            if line.startswith('#'):
                header += line
            else:
                sheader = header.split("\n")
                newinfo = sheader[-2] + "\tFORMAT\t"+ gto.header()
                sheader[-2] = newinfo
                header = "\n".join(sheader[:-1])
                outf.write(header + "\n")
                #print(line.rstrip())
                break
        for line in f:
            if line.startswith('#'):
                continue
            else:
                sline = line.rstrip().split()
                genotypes = gto.get_genotypes(sline[0], sline[7], ids_per_chrom)
                outf.write(line.rstrip() + "\tGT\t" + genotypes + "\n")
           



