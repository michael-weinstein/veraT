#!/usr/bin/env python3

class CheckArgs():  #class that checks arguments and ultimately returns a validated set of arguments to the main program
    
    def __init__(self):
        import argparse
        import os
        parser = argparse.ArgumentParser()
        parser.add_argument("-f", "--variantPickleFile", help = "Input pickle of fused, sorted variants.", required = True)
        parser.add_argument("-o", "--outputVariantFile", help = "Output file name", required = True)
        rawArgs = parser.parse_args()
        variantPickleFile = rawArgs.variantPickleFile
        if not os.path.isfile(variantPickleFile):
            raise FileNotFoundError("Unable to find input pickle file at %s" %variantPickleFile)
        self.variantPickleFile = variantPickleFile
        self.outputVariantFile = rawArgs.outputVariantFile

def makeOutputForOncotator(variantDict, delimiter = "\t"):
    import variantDataHandler
    variantList = list(variantDict.keys())
    variantDataHandler.sortVariantDataTuples(variantList)
    outputHeader = ["chr", "start", "end", "ref_allele", "alt_allele", "DNAcoverage", "Mut_DNA_Reads_Abs", "Mut_DNA_VAF", "NormalDNAcoverage"]
    outputLines = []
    outputLines.append(delimiter.join(outputHeader))
    for variant in variantList:
        if variantDict[variant].isIndel:
            continue
        outputLines.append(variantDict[variant].oncotatorInputLine(delimiter))
    outputLines = "\n".join(outputLines)
    return outputLines
    
if __name__ == '__main__':
    args = CheckArgs()
    import pickle
    variantPickleFile = open(args.variantPickleFile, 'rb')
    variantDict = pickle.load(variantPickleFile)
    variantPickleFile.close()
    fusedVariants = variantDict["fused"]["variants"]
    outputFile = open(args.outputVariantFile, 'w')
    outputFile.write(makeOutputForOncotator(fusedVariants))
    outputFile.close()
    quit()
    
