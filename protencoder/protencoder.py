from Bio import SeqIO
import numpy as np


class encoder(object):
    def __init__(self):
        self.seqDict = dict()
        self.filter = list()
        self.seqKeys = list()

    def read_fasta(self, seqPath):
        fasta_sequences = SeqIO.parse(seqPath, 'fasta')
        for fasta in fasta_sequences:
            if fasta.id in self.filter or self.filter == []:
                self.seqDict[fasta.id] = str(fasta.seq)

    def read_GO(self, filePath):
        if self.filter != []:
            self.seqKeys = self.filter
        keysF = open(filePath)
        for line in keysF.readlines():
            line = line.rstrip("\n").split("\t")
            if line[0] in self.filter or self.filter == []:
                if line[0] not in self.seqDict:
                    self.seqDict[line[0]] = {'F': [], 'P': [], 'C': []}
                self.seqDict[line[0]][line[2]].append(line[1])
                if self.filter == []:
                    self.seqKeys.append(line[0])
        keysF.close()

    def dump(self, outPrefix):
        seqArr = np.array(list(self.seqDict.values()))
        if self.filter == []:
            self.seqKeys = list(self.seqDict.keys())
        else:
            self.seqKeys = self.filter
        self.seqDict = dict()
        np.save(outPrefix + "_onehot", seqArr)
        del seqArr
        seqKeysF = open(outPrefix + "_keys.txt", 'w')
        for i in self.seqKeys:
            seqKeysF.write(i + "\n")
        seqKeysF.close()

    def dump_GO(self, GOclasses, outPrefix):
        ForderedGOA = []
        PorderedGOA = []
        CorderedGOA = []
        for seq in self.seqKeys:
            ForderedGOA.append(self.seqDict[seq]['F'])
            PorderedGOA.append(self.seqDict[seq]['P'])
            CorderedGOA.append(self.seqDict[seq]['C'])
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
        FGOAf = open(outPrefix + "_FGOclasses.txt", 'w')
        PGOAf = open(outPrefix + "_PGOclasses.txt", 'w')
        CGOAf = open(outPrefix + "_CGOclasses.txt", 'w')
        for GOA in GOclasses['F']:
            FGOAf.write(GOA + "\n")
        for GOA in GOclasses['P']:
            PGOAf.write(GOA + "\n")
        for GOA in GOclasses['C']:
            CGOAf.write(GOA + "\n")
        FGOAf.close()
        PGOAf.close()
        CGOAf.close()

    def load_filter(self, filter):
        seqKeysF = open(filter)
        for i in seqKeysF.readlines():
            self.filter.append(i.rstrip("\n"))
        seqKeysF.close()
