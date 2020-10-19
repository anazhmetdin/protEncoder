"""Console script for protencoder."""
import argparse
import sys
from protencoder.onehot import ONEencoder


def main():
    """Console script for protencoder."""
    parser = argparse.ArgumentParser()
    parser.add_argument('_', nargs='*')
    parser.add_argument("-d", "--seqPath", required=True,
                        help="file path of the protein sequences")
    parser.add_argument("-m", "--maxLen", default="2000",
                        help="maximum length of the protein to be encoded;\
                        default = 2000")
    parser.add_argument("-o", "--outPrefix", default="",
                        help="output files prefix")
    args = parser.parse_args()

    seqPath = args.seqPath
    maxLen = int(args.maxLen)
    outPrefix = args.outPrefix if args.outPrefix != "" else seqPath[:-6]

    oneHotencd = ONEencoder(maxLen)
    oneHotencd.handler.read_fasta(seqPath)
    oneHotencd.encode()
    oneHotencd.handler.dump(outPrefix)
    return 0


if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
