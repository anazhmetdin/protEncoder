from math import inf
from protencoder.protencoder import encoder


class ONEencoder():
    def __init__(self, maxLen):
        self.handler = encoder()
        self.maxLen = maxLen
        # Nine physicochemical properties for 20 amino acid types
        """
        H1, hydrophobicity
        H2, hydrophicility
        H3, hydrogen bond
        V, volumes of side chains
        P1, polarity
        P2, polarizability
        SASA, solvent-accessible surface area
        NCI, net charge index of side chains
        MASS, average mass of amino acid.
        """
        # This data is copied from: Li Z, Tang J, Guo F (2016)
        # DOI:10.1371/journal.pone.0147467
        self.aaDict = {'A': [0.62, -0.5, 2, 27.5, 8.1, 0.046, 1.181, 0.007187, 71.0788],
                       'C': [0.29, -1.0, 2.0, 44.6, 5.5, 0.128, 1.461, -0.03661, 103.1388],
                       'D': [-0.9, 3.0, 4.0, 40.0, 13.0, 0.105, 1.587, -0.02382, 115.0886],
                       'E': [-0.74, 3.0, 4.0, 62.0, 12.3, 0.151, 1.862, 0.006802, 129.1155],
                       'F': [1.19, -2.5, 2.0, 115.5, 5.2, 0.29, 2.228, 0.037552, 147.1766],
                       'G': [0.48, 0.0, 2.0, 0.0, 9.0, 0.0, 0.881, 0.179052, 57.0519],
                       'H': [-0.4, -0.5, 4.0, 79.0, 10.4, 0.23, 2.025, -0.01069, 137.1411],
                       'I': [1.38, -1.8, 2.0, 93.5, 5.2, 0.186, 1.81, 0.021631, 113.1594],
                       'K': [-1.5, 3.0, 2.0, 100.0, 11.3, 0.219, 2.258, 0.017708, 128.1741],
                       'L': [1.06, -1.8, 2.0, 93.5, 4.9, 0.186, 1.931, 0.051672, 113.1594],
                       'M': [0.64, -1.3, 2.0, 94.1, 5.7, 0.221, 2.034, 0.002683, 131.1986],
                       'N': [-0.78, 2.0, 4.0, 58.7, 11.6, 0.134, 1.655, 0.005392, 114.1039],
                       'P': [0.12, 0.0, 2.0, 41.9, 8.0, 0.131, 1.468, 0.239531, 97.1167],
                       'Q': [-0.85, 0.2, 4.0, 80.7, 10.5, 0.18, 1.932, 0.049211, 128.1307],
                       'R': [-2.53, 3.0, 4.0, 105.0, 10.5, 0.18, 1.932, 0.049211, 156.1875],
                       'S': [-0.18, 0.3, 4.0, 29.3, 9.2, 0.062, 1.298, 0.004627, 87.0782],
                       'T': [-0.05, -0.4, 4.0, 51.3, 8.6, 0.108, 1.525, 0.003352, 101.1051],
                       'V': [1.08, -1.5, 2.0, 71.5, 5.9, 0.14, 1.645, 0.057004, 99.1326],
                       'W': [0.81, -3.4, 3.0, 145.5, 5.4, 0.409, 2.663, 0.037977, 186.2132],
                       'Y': [0.26, -2.3, 3.0, 117.3, 6.2, 0.298, 2.368, 0.023599, 163.176]}
        self.Nprops = len(self.aaDict['A'])
        self.aaCount = len(self.aaDict)
        self.vecLen = self.Nprops + self.aaCount
        # creating an index for each AA
        self.aa = {b: [a] for a, b in enumerate(self.aaDict.keys())}
        # getting minimum and maximum value of each property
        MinMax = self.Nprops*[[inf, -inf]]
        for k in self.aaDict:
            for v in range(len(self.aaDict[k])):
                if self.aaDict[k][v] < MinMax[v][0]:
                    MinMax[v][0] = self.aaDict[k][v]
                if self.aaDict[k][v] > MinMax[v][1]:
                    MinMax[v][1] = self.aaDict[k][v]
        # normalize values of each property
        for k in self.aaDict:
            for v in range(len(self.aaDict[k])):
                self.aaDict[k][v] = self.aaDict[k][v] - MinMax[v][0]
                self.aaDict[k][v] /= (MinMax[v][1] - MinMax[v][0])
        # adding 'X', 'U', and 'O' for ambigous amino acids, Selenocysteine
        # and Pyrrolysine.
        self.aaDict['X'] = self.Nprops * [0.5]
        self.aa['X'] = list(range(len(self.aaDict)))
        self.aaDict['U'] = self.aaDict['C']
        self.aa['U'] = self.aa['C']
        self.aaDict['O'] = self.aaDict['K']
        self.aa['O'] = self.aa['K']
        # adding 'B' for 'N'/'D', and 'Z' for 'Q'/'E', and 'J' for 'L'/'I'
        summ = [sum(i) for i in zip(self.aaDict['N'], self.aaDict['D'])]
        self.aaDict['B'] = [x/2 for x in summ]
        self.aa['B'] = self.aa['N'] + self.aa['D']
        summ = [sum(i) for i in zip(self.aaDict['Q'], self.aaDict['E'])]
        self.aaDict['Z'] = [x/2 for x in summ]
        self.aa['Z'] = self.aa['Q'] + self.aa['E']
        summ = [sum(i) for i in zip(self.aaDict['L'], self.aaDict['I'])]
        self.aaDict['J'] = [x/2 for x in summ]
        self.aa['J'] = self.aa['L'] + self.aa['I']

    def encode(self):
        if self.handler.seqKeys != []:
            protList = list(self.handler.seqDict.keys())
            for prot in protList:
                if prot not in self.handler.seqDict:
                    self.handler.seqDict.pop(prot)
        for prot in self.handler.seqDict:
            encoded = []
            seq = self.handler.seqDict[prot]
            if len(seq) > self.maxLen:
                seq = seq[:self.maxLen]
            for ris in seq:
                zeros = len(self.aa)*[0]
                for pos in self.aa[ris]:
                    zeros[pos] = 1
                encoded += zeros
                encoded += self.aaDict[ris]
            if len(encoded) < self.maxLen*(self.vecLen):
                encoded += (self.maxLen*(29) - len(encoded))*[0]
            self.handler.seqDict[prot] = encoded

    def read(self, seqPath):
        self.handler.read_fasta(seqPath)
        if self.maxLen == -1:
            self.maxLen = len(sorted(self.handler.seqDict.values(),
                                     key=len)[-1])

    def dump(self, outPrefix):
        self.handler.dump(outPrefix)

    def load_filter(self, filter):
        self.handler.load_filter(filter)
