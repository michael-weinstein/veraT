#!/usr/bin/env python3

class CheckArgs():  #class that checks arguments and ultimately returns a validated set of arguments to the main program
    
    def __init__(self):
        import argparse
        import os
        parser = argparse.ArgumentParser()
        parser.add_argument("-c", "--hlaCall", help = "Athlates output for the case.  Can have multiple entries, all flagged as -c or --hlaCall.  Should be formatted with the molecule and then file, separated by a colon (A:/path/to/hlaCalls/A.typing.txt)", dest = "hlaCallFiles", action = "append", required = True)
        parser.add_argument("-v", "--variantsFile", help = "Pickle containing the variant calls and peptides for the case.", required = True)
        parser.add_argument("-o", "--output", help = "Output file name", required = True)
        rawArgs = parser.parse_args()
        hlaCallFiles = rawArgs.hlaCallFiles
        hlaCallFiles = [item.split(":") for item in hlaCallFiles]
        self.hlaCallFiles = {}
        hlaFileSet = set()
        for file in hlaCallFiles:
            if not len(file) == 2:
                raise RuntimeError("Each file should be passed as the molecule and file path separated by a colon. One of the passed values was:\n%s" %(":".join(file)))
            molecule, filePath = file
            molecule = molecule.upper()
            if not os.path.isfile(filePath):
                raise FileNotFoundError("Unable to find HLA call file at %s" %filepath)
            if not molecule in self.hlaCallFiles and not filePath in hlaFileSet:
                self.hlaCallFiles[molecule] = filePath
                hlaFileSet.add(filePath)
            elif molecule in self.hlaCallFiles:
                raise RuntimeError("HLA-%s call file was specified at least twice in arguments.\nFirst instance: %s\nSecond instance: %s" %(molecule, self.hlaCallFiles[molecule], filePath))
            elif filePath in hlaFileSet:
                raise RuntimeError("HLA call file %s was given at least twice in arguments. Each call file should only be listed once.  Please ensure that you did not list the same file for two different HLA molecules.")
        variantsFile = rawArgs.variantsFile
        if not os.path.isfile(variantsFile):
            raise FileNotFoundError("Unable to find variants pickle file at %s" %variantsFile)
        self.variantsFile = variantsFile
        self.output = rawArgs.output
        
def sortAllelesForCanonical(alleleSet):
    import operator
    alleleSet = [item.split(":") for item in alleleSet]
    alleleSet.sort(key = operator.itemgetter(1))
    alleleSet = [":".join(item) for item in alleleSet]
    return alleleSet
        
def readHLATypesFromAthlatesOutput(athlatesOutputFile):
    callFile = open(athlatesOutputFile, 'r')
    callData = callFile.read()
    callFile.close()
    calls = callData.split("Inferred Allelic Pairs")
    assert len(calls) == 2, "HLA call file %s appears to be malformed and is either missing an Inferred Allelic Pairs section or has multiples. Neither one should happen." %athlatesOutputFile
    calls = calls[1]
    calls = calls.split("\n")[2:]  #throwing out the first two lines from this.  The first should be the left-over hyphens from the title and the second should be a blank line for whitespace
    hlaCallLines = []
    for line in calls:
        line = line.strip()
        if line:
            hlaCallLines.append(line)
    hlaCallLines = [line.split() for line in hlaCallLines]
    predictedAlleleSet = set()
    predictedAlleleOrder = []
    for line in hlaCallLines:
        for allele in line[:-1]: #last column is not HLA type calls
            trimmedAllele = ":".join(allele.split(":")[:2])
            predictedAlleleSet.add(trimmedAllele)
            predictedAlleleOrder.append(trimmedAllele)
    predictedAlleleOrder = sortAllelesForCanonical(predictedAlleleOrder)
    if len(predictedAlleleSet) < 3:
        return list(predictedAlleleSet)
    else:
        topAlleles = set()
        for allele in predictedAlleleOrder:
            topAlleles.add(allele)
            if len(topAlleles) == 2:
                return list(topAlleles)
            
def main():
    import pickle
    args = CheckArgs()
    variantPickleFile = open(args.variantsFile, 'rb')
    variantPickle = pickle.load(variantPickleFile)
    variantPickleFile.close()
    if "hla" in variantPickle:
        raise RuntimeError("HLA calls appear to already be present in this variant pickle.")
    hlaMolecules = list(args.hlaCallFiles.keys())
    hlaCalls = {}
    for molecule in hlaMolecules:
        hlaCalls[molecule] = readHLATypesFromAthlatesOutput(args.hlaCallFiles[molecule])
    variantPickle["hla"] = hlaCalls.copy()
    outputFile = open(args.output, 'wb')
    pickle.dump(variantPickle, outputFile)
    outputFile.close()
    
if __name__ == '__main__':
    main()
    