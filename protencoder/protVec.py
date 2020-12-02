import biovec
import random
from protencoder.protencoder import encoder


class protvec():
    def __init__(self, model_path="swissprot-reviewed-protvec.model",
                 flatten=False):
        self.handler = encoder()
        self.pv = biovec.models.load_protvec(model_path)
        self.flatten = flatten

    def encode(self):
        encoded = []
        Xs = ['R', 'K', 'D', 'Q', 'N', 'E', 'H', 'S', 'T',
              'P', 'Y', 'C', 'G', 'A', 'M', 'W', 'L', 'V',
              'F', 'I']
        Bs = ['N', 'D']
        Zs = ['Q', 'E']
        Js = ['L', 'I']
        for prot in self.handler.seqDict:
            seq = self.handler.seqDict[prot]
            while ('X' in seq):
                aa = Xs[random.randint(0, len(Xs)-1)]
                seq.replace("X", aa, 1)
            while ('B' in seq):
                aa = Bs[random.randint(0, len(Bs)-1)]
                seq.replace("B", aa, 1)
            while ('Z' in seq):
                aa = Zs[random.randint(0, len(Zs)-1)]
                seq.replace("Z", aa, 1)
            while ('J' in seq):
                aa = Js[random.randint(0, len(Js)-1)]
                seq.replace("J", aa, 1)
            if 'U' in seq:
                seq.replace("U", 'C')
            if 'O' in seq:
                seq.replace("O", "k")
            encoded = self.pv.to_vecs(seq)
            if self.flatten:
                encoded = encoded.flatten()
            self.handler.seqDict[prot] = encoded

    def read(self, seqPath):
        self.handler.read_fasta(seqPath)

    def dump(self, outPrefix):
        self.handler.dump(outPrefix, "protVec")
