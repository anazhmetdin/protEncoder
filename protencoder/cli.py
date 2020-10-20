"""Console script for protencoder."""
import argparse
import sys
from protencoder.onehot import ONEencoder
from protencoder.GOencoder import GOencoder


def main():
    """Console script for protencoder."""
    parser = argparse.ArgumentParser()
    parser.add_argument('_', nargs='*')
    filePathParser = parser.add_mutually_exclusive_group(required=True)
    filePathParser.add_argument("-d", "--seqPath",
                                help="file path of the protein sequences")
    filePathParser.add_argument("-g", "--GOfile",
                                help="file path of the GOA of proteins")
    parser.add_argument("-f", "--filter", default="",
                        help="file path of protein keys to be used for\
                        encoding filteration")
    parser.add_argument("-m", "--maxLen", default="2000",
                        help="maximum length of the protein to be encoded;\
                        if -1 m = max protein length, default = 2000")
    parser.add_argument("-o", "--outPrefix", default="",
                        help="output files prefix")
    args = parser.parse_args()

    seqPath = args.seqPath
    GOfile = args.GOfile
    filter = args.filter
    maxLen = int(args.maxLen)

    if not (seqPath is None):
        outPrefix = args.outPrefix if args.outPrefix != "" else seqPath[:-6]
        oneHotencd = ONEencoder(maxLen)
        if filter != "":
            oneHotencd.load_filter(filter)
        oneHotencd.read(seqPath)
        oneHotencd.encode()
        oneHotencd.dump(outPrefix)

    elif not (GOfile is None):
        outPrefix = args.outPrefix if args.outPrefix != "" else GOfile[:-4]
        oneHotencd = GOencoder()
        if filter != "":
            oneHotencd.load_filter(filter)
        oneHotencd.read(GOfile)
        oneHotencd.encode()
        oneHotencd.dump(outPrefix)
    return 0


if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
