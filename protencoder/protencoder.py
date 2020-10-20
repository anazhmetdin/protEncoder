from Bio import SeqIO
import numpy as np


class encoder(object):
    def __init__(self):
        self.seqDict = dict()
        self.seqKeys = list()

    def read_fasta(self, seqPath):
        fasta_sequences = SeqIO.parse(seqPath, 'fasta')
        for fasta in fasta_sequences:
            if fasta.id in self.seqKeys or self.seqKeys == []:
                self.seqDict[fasta.id] = str(fasta.seq)

    def read_GO(self, filePath):
        keysF = open(filePath)
        for line in keysF.readlines():
            line = line.rstrip("\n").split("\t")
            if line[0] in self.seqKeys or self.seqKeys == []:
                if line[0] not in self.seqDict:
                    self.seqDict[line[0]] = {'F': [], 'P': [], 'C': []}
                self.seqDict[line[0]][line[2]].append(line[1])

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

    def dump_GO(self, outPrefix):
        ForderedGOA = []
        PorderedGOA = []
        CorderedGOA = []
        for seq in self.seqKeys:
            ForderedGOA.append([self.seqDict[seq]['F']])
            PorderedGOA.append([self.seqDict[seq]['P']])
            CorderedGOA.append([self.seqDict[seq]['C']])
            self.seqDict.pop(seq)
        ForderedGOA = np.array(ForderedGOA)
        PorderedGOA = np.array(PorderedGOA)
        CorderedGOA = np.array(CorderedGOA)
        np.save(outPrefix + "_FGOA", ForderedGOA)
        del ForderedGOA
        np.save(outPrefix + "_PGOA", PorderedGOA)
        del PorderedGOA
        np.save(outPrefix + "_CGOA", CorderedGOA)
        del CorderedGOA

    def load(self, filter):
        seqKeysF = open(filter)
        for i in seqKeysF.readlines():
            self.seqKeys.append(i.rstrip("\n"))
