from Bio import SeqIO
import numpy as np


class encoder(object):
    def __init__(self):
        self.seqDict = dict()
        self.seqKeys = list()

    def read_fasta(self, seqPath):
        fasta_sequences = SeqIO.parse(seqPath, 'fasta')
        for fasta in fasta_sequences:
            self.seqDict[fasta.id] = str(fasta.seq)

    def dump(self, outPrefix):
        seqArr = np.array(list(self.seqDict.values()))
        seqKeys = list(self.seqDict.keys())
        self.seqDict = dict()
        np.save(outPrefix + "_onehot", seqArr)
        del seqArr
        seqKeysF = open(outPrefix + "_keys.txt", 'w')
        for i in seqKeys:
            seqKeysF.write(i + "\n")
        seqKeysF.close()

    def load(self, inPreFix):
        seqKeysF = open(inPreFix + "_keys.txt")
        for i in seqKeysF.readline():
            self.seqKeys.append(i.rstrip("\n"))
