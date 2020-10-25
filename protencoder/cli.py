"""Console script for protencoder."""
import argparse
import sys
import glob
from Bio import SeqIO
from collections import Counter
from protencoder.onehot import ONEencoder
from protencoder.GOencoder import GOencoder


def devide(seqPath, chopSize, outPrefix):
    fasta_sequences = SeqIO.parse(seqPath, 'fasta')
    x = 0
    newFasta = list()
    part = 1
    for seq in fasta_sequences:
        newFasta.append(seq)
        x += 1
        if x == chopSize:
            SeqIO.write(newFasta, outPrefix+"_part"+str(part)+".fasta",
                        'fasta')
            newFasta = list()
            x = 0
            part += 1
    SeqIO.write(newFasta, outPrefix+"_part"+str(part)+".fasta",
                'fasta')


def create_filter(collection, numFreqGO, outPrefix):
    seqDict = {'F': list(), 'P': list(), 'C': list()}
    GOAf = open(collection)
    for line in GOAf.readlines():
        line = line.rstrip("\n").split("\t")
        seqDict[line[2]].append(line[1])
    GOAf.close()
    ofile = outPrefix+"_filter_"+str(numFreqGO)
    filterF = open(ofile, 'w')
    counter = Counter(seqDict['F'])
    most_occur = [x[0] for x in counter.most_common(numFreqGO)]
    for GOA in most_occur:
        filterF.write(GOA+"\t"+'F'+"\n")

    counter = Counter(seqDict['P'])
    most_occur = [x[0] for x in counter.most_common(numFreqGO)]
    for GOA in most_occur:
        filterF.write(GOA+"\t"+'P'+"\n")

    counter = Counter(seqDict['C'])
    most_occur = [x[0] for x in counter.most_common(numFreqGO)]
    for GOA in most_occur:
        filterF.write(GOA+"\t"+'C'+"\n")
    filterF.close()
    return ofile


def main():
    """Console script for protencoder."""
    parser = argparse.ArgumentParser()
    parser.add_argument('_', nargs='*')
    filePathParser = parser.add_mutually_exclusive_group(required=True)
    filePathParser.add_argument("-d", "--seqPath",
                                help="file path of the protein sequences")
    filePathParser.add_argument("-g", "--GOfile",
                                help="file path of the GOA of proteins")
    parser.add_argument("-f", "--Protfilter", default="",
                        help="file path of protein keys to be used for\
                        encoding filteration. in case of multiple filters\
                        enter the consensus path followed by *")
    parser.add_argument("-m", "--maxLen", default="2000",
                        help="maximum length of the protein to be encoded;\
                        if -1 m = max protein length, default = 2000")
    parser.add_argument("-p", "--GOpartioned", default="0",
                        help="if annotated keys are divided;\
                        0(default) if False, 1 if True")
    parser.add_argument("-o", "--outPrefix", default="",
                        help="output files prefix")
    parser.add_argument("-s", "--chopSize", default="-1",
                        help="number of sequences to be loaded\
                        into memory per time")
    GOfilterParser = parser.add_mutually_exclusive_group()
    GOfilterParser.add_argument("-F", "--GOfilter", default="",
                                help="file path of GO keys to be used for\
                                encoding filteration")
    GOfilterParser.add_argument("-c", "--collection", default="",
                                help="file path of collected\
                                protein annotations")
    parser.add_argument("-n", "--numFreqGO", default="500", help="number of\
                        the most frequent GO annotions to be considered")
    args = parser.parse_args()

    seqPath = args.seqPath
    GOfile = args.GOfile
    Protfilter = args.Protfilter
    GOfilter = args.GOfilter
    maxLen = int(args.maxLen)
    chopSize = int(args.chopSize)
    GOpartioned = bool(int(args.GOpartioned))
    collection = args.collection
    numFreqGO = int(args.numFreqGO)
    outPrefix = args.outPrefix

    if not (seqPath is None):
        outPrefix = args.outPrefix if args.outPrefix != "" else seqPath[:-6]
        if collection != "":
            create_filter(collection, numFreqGO, outPrefix)
        oneHotencd = ONEencoder(maxLen)
        if Protfilter != "":
            oneHotencd.load_filter(Protfilter)
        if chopSize == -1:
            oneHotencd.read(seqPath)
            oneHotencd.encode()
            oneHotencd.dump(outPrefix)
        else:
            devide(seqPath, chopSize, outPrefix)
            for name in glob.glob(outPrefix+"_part*.fasta"):
                oneHotencd.read(name)
                oneHotencd.encode()
                num = name[name.find("part") + 4:-6]
                oneHotencd.dump(outPrefix+"_part"+num)

    elif not (GOfile is None):
        outPrefix = args.outPrefix if args.outPrefix != "" else GOfile[:-4]
        if collection != "":
            GOfilter = create_filter(collection, numFreqGO, outPrefix)
        if not GOpartioned:
            oneHotencd = GOencoder()
            if Protfilter != "":
                oneHotencd.load_filter(Protfilter, GOfilter)
            oneHotencd.read(GOfile)
            oneHotencd.encode()
            oneHotencd.dump(outPrefix)
        else:
            for name in glob.glob(Protfilter):
                oneHotencd = GOencoder()
                oneHotencd.load_filter(name, GOfilter)
                oneHotencd.read(GOfile)
                oneHotencd.encode()
                num = name[name.find("part") + 4:-6]
                oneHotencd.dump(outPrefix+"_part"+num)
    return 0


if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
