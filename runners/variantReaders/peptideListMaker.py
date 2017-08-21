#!/usr/bin/env python3

class CheckArgs():  #class that checks arguments and ultimately returns a validated set of arguments to the main program
    
    def __init__(self):
        import argparse
        import os
        parser = argparse.ArgumentParser()
        parser.add_argument("-v", "--variantsFile", help = "Pickle containing the variant calls and peptides for the case.", required = True)
        parser.add_argument("-o", "--output", help = "Output base file name", required = True)
        rawArgs = parser.parse_args()
        variantsFile = rawArgs.variantsFile
        if not os.path.isfile(variantsFile):
            raise FileNotFoundError("Unable to find variants pickle file at %s" %variantsFile)
        self.variantsFile = variantsFile
        self.output = rawArgs.output
        
def makePeptideDict(fusedVariantDict):
    peptideDict = {}
    for variantHash in list(fusedVariantDict['neoepitopes'].keys()):
        for peptide in fusedVariantDict['neoepitopes'][variantHash]:
            if not len(peptide) in peptideDict:
                peptideDict[len(peptide)] = set()
            peptideDict[len(peptide)].add(peptide)
        for peptide in fusedVariantDict['wildtype'][variantHash]:
            peptideDict[len(peptide)].add(peptide)
    return (peptideDict, list(peptideDict.keys()))
             
def makePeptideList(args):
    import pickle
    variantPickleFile = open(args.variantsFile, 'rb')
    variantPickle = pickle.load(variantPickleFile)
    variantPickleFile.close()
    peptideDict, determinedLengths = makePeptideDict(variantPickle["fused"])  #length list appears to not load properly on the cluster.  Not sure why, same file loads properly locally with the correct length list.  Made sure variantDataHandler version were identical and MD5 values match between pickles.  Potentially a bug somewhere?
    peptideLengths = sorted(list(peptideDict.keys()))
    if not peptideLengths:
        peptideLengths = sorted(determinedLengths)
    variantPickle["peptideLengths"] = peptideLengths
    peptideFileDict = {}
    for length in peptideLengths:
        outputFileName = args.output + ".%s.pep" %length
        peptideFileDict[length] = outputFileName
        outputFile = open(outputFileName, 'w')
        outputFile.write("\n".join(list(peptideDict[length])))
        outputFile.close()
    variantPickleOutput = open(args.variantsFile, 'wb')
    pickle.dump(variantPickle, variantPickleOutput)
    variantPickleOutput.close()
    return peptideFileDict
    
if __name__ == '__main__':
    args = CheckArgs()
    makePeptideList(args)
    