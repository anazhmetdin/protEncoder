from protencoder.protencoder import encoder
import random
from itertools import product
from copy import deepcopy


class protKmers():
    def __init__(self, k):
        self.handler = encoder()
        self.k = k
        self.aa = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N',
                   'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
        self.kmers = list(product(self.aa, repeat=k))
        self.encodedTemp = [0 for i in self.kmers]

    def encode(self):
        encoded = []
        for prot in self.handler.seqDict:
            encoded = deepcopy(self.encodedTemp)
            seq = self.handler.seqDict[prot]
            while ('X' in seq):
                aa = Xs[random.randint(0, len(Xs)-1)]
                seq = seq.replace("X", aa, 1)
            while ('B' in seq):
                aa = Bs[random.randint(0, len(Bs)-1)]
                seq = seq.replace("B", aa, 1)
            while ('Z' in seq):
                aa = Zs[random.randint(0, len(Zs)-1)]
                seq = seq.replace("Z", aa, 1)
            while ('J' in seq):
                aa = Js[random.randint(0, len(Js)-1)]
                seq = seq.replace("J", aa, 1)
            if 'U' in seq:
                seq = seq.replace("U", 'C')
            if 'O' in seq:
                seq = seq.replace("O", "K")
            for i in range(len(self.kmers)):
                kmer = "".join(list(self.kmers[i]))
                freq = seq.count(kmer)
                encoded[i] = freq
            self.handler.seqDict[prot] = encoded

    def read(self, seqPath):
        self.handler.read_fasta(seqPath)

    def dump(self, outPrefix):
        self.handler.dump(outPrefix, "kmerHz")
