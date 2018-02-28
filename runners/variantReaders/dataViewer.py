#! /usr/bin/python3

fileName = "LB2907_DNA.acceptedSomatics.mPileupRNASupport.pkl"
if not fileName:
    raise RuntimeError("File name has to be set")

import pickle
file = open(fileName, 'rb')
data = pickle.load(file)
file.close()
variantsOfInterest = {}
for variant in data:
    if data[variant]["RNASupport"].totalDepth:
        variantsOfInterest[variant] = data[variant]
print(fileName)