#!/usr/bin/env python3

class CheckArgs(object):
    
    def __init__(self):
        import os
        import argparse
        parser = argparse.ArgumentParser()
        parser.add_argument("-i", "--inputFasta", help = "Input fasta to convert", required = True)
        parser.add_argument("-o", "--outputPickle", help = "Where to output the pickle containing the peptide data", required = True)
        parser.add_argument("-9", "--clobber", help = "Overwrite existing file with new data", action = 'store_true')
        rawArgs = parser.parse_args()
        inputFasta = rawArgs.inputFasta
        if not os.path.isfile(inputFasta):
            raise FileNotFoundError("Unable to find fasta file at %s" %inputFasta)
        self.inputFasta = inputFasta
        outputPickle = rawArgs.outputPickle
        if not rawArgs.clobber and os.path.isfile(outputPickle):
            raise FileExistsError("Output file %s already exists. Please choose a different file name, move or delete the original or rerun this program in clobber mode (adding a -9 to the arguments)")
        self.outputPickle = outputPickle

class ProteinFastaLine(object):
    
    def __init__(self, line):
        if not line:
            raise RuntimeError("Empty line passed as protein fasta line.")
        self.line = line.strip()
        lineArray = self.line.split("\n")
        lineData = lineArray[0]
        self.peptideSequence = "".join(lineArray[1:])
        infoArray = lineData.split()
        self.ensp = infoArray[0]
        del infoArray[0]
        infoArray = [item.split(":") for item in infoArray]
        infoTable = {}
        for item in infoArray:
            infoTable[item[0]] = ":".join(item[1:])
        self.transcript = infoTable["transcript"]
        self.gene = infoTable["gene"]

def deleteEmptyListEntries(listToProcess):  #this will do the deletions in place
    for i in range(len(listToProcess) - 1, -1, -1):
        if not listToProcess[i]:
            del listToProcess[i]

def analyzeProteinFasta(fileName):
    if fileName.endswith(".gz"):
        import gzip
        fastaFile = gzip.open(fileName, 'rt')
    else:
        fastaFile = open(fileName, 'r')
    gene2Transcript = {}
    transcript2Protein = {}
    protein2Sequence = {}
    fasta = fastaFile.read()
    fastaFile.close()
    fasta = fasta.split(">")
    deleteEmptyListEntries(fasta)
    for entry in fasta:
        data = ProteinFastaLine(entry)
        if not data.gene in gene2Transcript:
            gene2Transcript[data.gene] = []
        gene2Transcript[data.gene].append(data.transcript)
        transcript2Protein[data.transcript] = data.ensp
        protein2Sequence[data.ensp] = data.peptideSequence
    return {'gene2Transcript' : gene2Transcript,
            'transcript2Protein' : transcript2Protein,
            'protein2Sequence' : protein2Sequence}

class PeptideSequenceGetter(object):
    
    def __init__(self, peptideSequencePickle):
        import pickle
        file = open(peptideSequencePickle, 'rb')
        self.peptideData = pickle.load(file)
        file.close()
        
    def getPeptide(self, identifier):
        identifier = identifier.upper().split(".")[0]  #make it uppercase and get rid of any dot values after the id number (they cause mismatches)
        if identifier.startswith("ENSP"):
            return self.getPeptideFromENSP(identifier)
        elif identifier.startswith("ENST"):
            return self.getPeptideFromENST(identifier)
        elif identifier.startswith("ENSG"):
            return self.getPeptideFromENSG(identifier)
        else:
            raise RuntimeError("Peptide must be found from either a peptide, transcript, or gene identifier (ENSP/T/G number)")
        
    def getPeptideFromENSP(self, identifier):
        return self.peptideData["protein2Sequence"][identifier]
        
    def getPeptideFromENST(self, identifier):
        identifier = identifier.split(".")[0]
        return self.getPeptideFromENSP(self.peptideData["transcript2Protein"][identifier])
        
    def getPeptideFromENSG(self, identifier):
        transcriptList = self.peptideData["gene2Transcript"][identifier]
        peptideList = []
        for transcript in transcriptList:
            peptideList.append(self.getPeptideFromENST(transcript), transcript)
        return peptideList
    
def createPeptidePickle(inputFasta, outputPickle):
    import pickle
    peptidePickle = analyzeProteinFasta(inputFasta)
    outputFile = open(outputPickle, 'wb')
    pickle.dump(peptidePickle, outputFile)
    outputFile.close()
    
# getter = PeptideSequenceGetter("homoSapiens.pkl")
# test = getter.getPeptide("ENST00000418703")  #using TRPV4 here because I know where a few residues should be
# r = test[592:595]
# p = test[797:800]
# print("something")
    
if __name__ == '__main__':
    import datetime
    start = datetime.datetime.now()
    args = CheckArgs()
    createPeptidePickle(args.inputFasta, args.outputPickle)
    print("Pickle created in %s" %(datetime.datetime.now() - start))
    quit()