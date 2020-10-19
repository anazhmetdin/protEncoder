from Bio import SeqIO
import numpy as np


class encoder(object):
    def __init__(self):
        self.seqDict = dict()

    def read_fasta(self, seqPath):
        fasta_sequences = SeqIO.parse(seqPath, 'fasta')
        for fasta in fasta_sequences:
            self.seqDict[fasta.id] = str(fasta.seq)

    def dump(self, outPrefix):
        seqArr = np.array(list(self.seqDict.values()))
        ordered = list(self.seqDict.keys())
        self.seqDict = dict()
        np.save(outPrefix + "_onehot", seqArr)
        del seqArr
        orderedF = open(outPrefix + "_ordered", 'w')
        for i in ordered:
            orderedF.write(i + "\n")
