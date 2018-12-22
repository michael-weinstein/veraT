#! /usr/bin/python3

def openVariantsPickle(fileName):
    if not fileName:
        raise RuntimeError("File name has to be set")
    import pickle
    file = open(fileName, 'rb')
    data = pickle.load(file)
    file.close()
    return data

if __name__ == "__main__":
    fileName = "DF12GD_tumorgDNA.acceptedSomatics.mPileupRNASupport.fused.oncotatorAnnotation.peptides.hlaCalls.pkl"
    data = openVariantsPickle(fileName)
    print("Data headers: %s" %data.keys())