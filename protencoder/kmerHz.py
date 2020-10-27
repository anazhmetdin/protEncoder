from protencoder.protencoder import encoder
from itertools import combinations_with_replacement
from copy import deepcopy


class protKmers():
    def __init__(self, k):
        self.handler = encoder()
        self.k = k
        self.aa = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N',
                   'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
        self.kmers = list(combinations_with_replacement(self.aa, k))
        self.encodedTemp = [0 for i in self.kmers]

    def encode(self):
        encoded = []
        for prot in self.handler.seqDict:
            encoded = deepcopy(self.encodedTemp)
            seq = self.handler.seqDict[prot]
            for i in range(len(self.kmers)):
                kmer = "".join(list(self.kmers[i]))
                freq = seq.count(kmer)
                encoded[i] = freq
            self.handler.seqDict[prot] = encoded
        print(encoded)

    def read(self, seqPath):
        self.handler.read_fasta(seqPath)
