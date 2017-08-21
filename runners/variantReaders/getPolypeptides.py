#!/usr/bin/env python3

class CheckArgs(object):
    
    def __init__(self):
        import os
        import argparse
        parser = argparse.ArgumentParser()
        parser.add_argument("-v", "--variantPickle", help = "Pickle containing variants", required = True)
        parser.add_argument("-p", "--proteinSequencePickle", help = "Pickle containing protein sequence data and ensembl id conversions.  Make this by running the protein sequence finder module as an executable")
        parser.add_argument("-o", "--outputPickle", help = "Where to output the pickle containing the peptide lists for each variant", required = True)
        parser.add_argument("-s", "--sizeList", help = "Comma-separated list of desired peptide sizes", default = "8,9,10")
        rawArgs = parser.parse_args()
        variantPickle = rawArgs.variantPickle
        if not os.path.isfile(variantPickle):
            raise FileNotFoundError("Unable to find fasta file at %s" %variantPickle)
        self.variantPickle = variantPickle
        proteinSequencePickle = rawArgs.proteinSequencePickle
        if not os.path.isfile(proteinSequencePickle):
            raise FileNotFoundError("Unable to find fasta file at %s" %proteinSequencePickle)
        self.proteinSequencePickle = proteinSequencePickle
        self.outputPickle = rawArgs.outputPickle
        sizeList = rawArgs.sizeList
        sizeList = sizeList.split(",")
        self.sizeList = [int(item) for item in sizeList]

def matchingReferences(changes, proteinSequence):
    for change in changes:
        ref, pos, alt = change
        pos = pos - 1  #string is indexed to zero and variant position is indexed to 1.  This fixes that.
        if not ref == proteinSequence[pos]:
            return False
    return True

class TargetedWindow(object):
    
    def __init__(self, size, targets, sequence, targetOffset = -1, requireAllTargetsInWindow = False):  #target offset will be the difference in indexing between targets (which may be indexed to 1) and python strings (which may be indexed to 0).  Example: P19S is indexed to 1 (first amino acid = 1), while the first character of the string would be zero.  Offset will be -1.
        self.size = size
        if type(targets) == int:
            targets = [targets]
        self.targets = targets
        self.requireAllTargetsInWindow = requireAllTargetsInWindow
        self.targets = [target + targetOffset for target in targets]
        self.sequence = sequence
        self.end = self.size
        self.start = 0
        self.done = False
        
    def reset(self):
        self.end = self.size
        self.start = 0
        self.done = False
        
    def advance(self, positions = 1):
        self.start += positions
        self.end += positions
        self.done = self.end >= len(self.sequence)
        
    def getNextWindow(self):
        while self.start < 0: #this should be impossible
            self.advance()
        while not self.done and not self.onTarget():
            self.advance()
        if self.done:
            return False
        else:
            windowSeq =  self.sequence[self.start:self.end]
            self.advance()
            return windowSeq
            
    def onTarget(self):
        targetSet = set(self.targets)
        windowSet = set(range(self.start, self.end))
        onTargetSet = targetSet.intersection(windowSet)
        if self.requireAllTargetsInWindow:
            return onTargetSet == targetSet
        else:
            return bool(onTargetSet)
        
    def getAllWindows(self):
        self.reset()
        windows = []
        while not self.done:
            window = self.getNextWindow()
            if window:
                windows.append(window)
        return windows
        
def makeMutantProtein(proteinSequence, changes, targetOffset = -1):
    proteinSequence = list(proteinSequence)
    for change in changes:
        proteinSequence[change[1] + targetOffset] = change[2]
    return "".join(proteinSequence)

def unmatchedSequenceInfo(proteinSequence, changes, targetOffset = -1):
    returnData = []
    for change in changes:
        returnData.append("Expected change: %s%s%s" %(change[0], change[1], change[2]))
        returnData.append("Found reference: %s%s" %(proteinSequence[change[1] + targetOffset], change[1]))
    return "\n".join(returnData)
    
def getPeptides(proteinSequencePickle, proteinChangeDict, peptideSizes = [8,9,10], unmatchedSequenceHandling = "skip"):
    peptideListDict = {"wildtype":{}, "mutant":{}}
    import proteinSequenceFinder
    proteinGetter = proteinSequenceFinder.PeptideSequenceGetter(proteinSequencePickle)
    for variant in list(proteinChangeDict.keys()):
        currentPeptideList = []
        changeData = proteinChangeDict[variant]
        enst = changeData.enst
        changes = changeData.changes
        proteinSequence = proteinGetter.getPeptide(enst)
        mutationSites = [change[1] for change in changes]
        if not matchingReferences (changes, proteinSequence):
            if unmatchedSequenceHandling == "skip":
                failedMatchData = unmatchedSequenceInfo(proteinSequence, changes)
                print("Unmatching sequences in %s: \n%s\n%s" %(changeData.gene, failedMatchData, proteinSequence))
                peptideListDict[variant] = currentPeptideList
                continue
            elif unmatchedSequenceHandling == "fail":
                failedMatchData = unmatchedSequenceInfo(proteinSequence, changes)
                raise RuntimeError("Unmatching sequences in %s: \n%s\n%s" %(changeData.gene, failedMatchData, proteinSequence))
            elif unmatchedSequenceHandling == "force":
                pass
            else:
                raise ValueError("Found an unmatching reference amino acid and an invalid handling option was given. Option: %s" %unmatchedSequenceHandling)
                mutantProtein = makeMutantProtein(proteinSequence, changes)
        for peptideSize in peptideSizes:
            targetWindow = TargetedWindow(peptideSize, mutationSites, proteinSequence)
            currentPeptideList += targetWindow.getAllWindows()
        peptideListDict["wildtype"][variant] = currentPeptideList.copy()
        currentPeptideList = []
        mutantProtein = makeMutantProtein(proteinSequence, changes)
        #print(proteinSequence[mutationSites[0] - 1])
        #print(mutantProtein[mutationSites[0] - 1])
        for peptideSize in peptideSizes:
            targetWindow = TargetedWindow(peptideSize, mutationSites, mutantProtein)
            currentPeptideList += targetWindow.getAllWindows()
        peptideListDict["mutant"][variant] = currentPeptideList
    return peptideListDict

def getVariantPickle(variantPickleFileName):
    import pickle
    variantPickleFile = open(variantPickleFileName, 'rb')
    variantPickle = pickle.load(variantPickleFile)
    variantPickleFile.close()
    return variantPickle

def main():
    args = CheckArgs()
    import pickle
    variantPickle = getVariantPickle(args.variantPickle)
    peptideChanges = variantPickle["fused"]["proteinChanges"]
    peptideListDict = getPeptides(args.proteinSequencePickle, peptideChanges, args.sizeList)
    variantPickle["fused"]["wildtype"] = peptideListDict["wildtype"]
    variantPickle["fused"]["neoepitopes"] = peptideListDict["mutant"]
    outputFile = open(args.outputPickle, 'wb')
    pickle.dump(variantPickle, outputFile)
    outputFile.close()
    quit()
    
if __name__ == '__main__':
    main()