import biovec
from protencoder.protencoder import encoder


class protvec():
    def __init__(self, model_path="swissprot-reviewed-protvec.model",
                 flatten=False):
        self.handler = encoder()
        self.pv = biovec.models.load_protvec(model_path)
        self.flatten = flatten

    def encode(self):
        encoded = []
        for prot in self.handler.seqDict:
            seq = self.handler.seqDict[prot]
            encoded = self.pv.to_vecs(seq)
            if self.flatten:
                encoded = encoded.flatten()
            self.handler.seqDict[prot] = encoded

    def read(self, seqPath):
        self.handler.read_fasta(seqPath)

    def dump(self, outPrefix):
        self.handler.dump(outPrefix, "kmerHz")
