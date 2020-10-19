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
    args = parser.parse_args()

    seqPath = args.seqPath

    oneHotencd = ONEencoder()
    oneHotencd.handler.read_fasta(seqPath)
    oneHotencd.encode()
    oneHotencd.handler.dump(seqPath[:-4])
    return 0


if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
