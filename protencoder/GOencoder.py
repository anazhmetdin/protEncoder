from protencoder.protencoder import encoder


class GOencoder(encoder):
    def __init__(self):
        self.handler = encoder()
        self.GOclasses = {'F': set(), 'P': set(), 'C': set()}

    def encode(self):
        for prot in self.handler.seqDict:
            for GOclass in self.handler.seqDict[prot]:
                for GOA in self.handler.seqDict[prot][GOclass]:
                    self.GOclasses[GOclass].add(GOA)

        for GOclass in self.GOclasses:
            self.GOclasses[GOclass] = list(self.GOclasses[GOclass])
            self.GOclasses[GOclass] = {b: a for a,
                                       b in enumerate(self.GOclasses[GOclass])}
        for seq in self.handler.seqDict:
            for GOclass in self.handler.seqDict[seq]:
                encoded = len(self.GOclasses[GOclass]) * [0]
                for GOA in self.handler.seqDict[seq][GOclass]:
                    encoded[self.GOclasses[GOclass][GOA]] = 1
                self.handler.seqDict[seq][GOclass] = encoded

    def read(self, GOfile):
        self.handler.read_GO(GOfile)

    def dump(self, outPrefix):
        self.handler.dump_GO(self.GOclasses, outPrefix)

    def load_filter(self, filter):
        self.handler.load_filter(filter)
