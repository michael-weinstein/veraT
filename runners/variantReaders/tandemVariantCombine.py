#!/usr/bin/env python3

class CheckArgs():  #class that checks arguments and ultimately returns a validated set of arguments to the main program
    
    def __init__(self):
        import argparse
        import os
        parser = argparse.ArgumentParser()
        parser.add_argument("-f", "--variantPickleFile", help = "Input pickle of unfused, sorted variants.", required = True)
        parser.add_argument("-o", "--outputVariantFile", help = "Output file name", required = True)
        parser.add_argument("-m", "--maxFusionLength", help = "Maximum number of bases to fuse into a single variant (value of 1 will give max fusions of length 2, 3 will give 2, etc.  Value of 0 will set no maximum length.)", type = int, default = 0)
        parser.add_argument("-p", "--maxDifferencePercent" , help = "Maximum percent difference in variant read percentage between one locus and the next before it is rejected as a tandem. Enter values as a whole number (10 for 10 percent, 100 to not check at all).", type = int, default = 10)
        rawArgs = parser.parse_args()
        variantPickleFile = rawArgs.variantPickleFile
        if not os.path.isfile(variantPickleFile):
            raise FileNotFoundError("Unable to find input pickle file at %s" %variantPickleFile)
        self.variantPickleFile = variantPickleFile
        self.outputVariantFile = rawArgs.outputVariantFile
        maxFusionLength = rawArgs.maxFusionLength
        if maxFusionLength < 0:
            raise RuntimeError("Max fusion length must be a value of zero or greater. %s was given." %maxFusionLength)
        self.maxFusionLength = maxFusionLength
        maxDifferencePercent = rawArgs.maxDifferencePercent
        if maxDifferencePercent < 0:
            raise RuntimeError("Max difference percent must be a value of zero or greater. %s was given." %maxDifferencePercent)
        self.maxDifferencePercent = maxDifferencePercent

def makeCombinedAndRNAEntries(variantDict):
    keys = list(variantDict.keys())
    variantDict["combined"] = {}
    if "RNASupport" in variantDict[keys[0]]:
        usingRNA = True
    else:
        usingRNA = False
    if usingRNA:
        variantDict["RNASupport"] = {}
    for variantHash in keys:
        variantDict["combined"][variantHash] = variantDict[variantHash]["combined"]
        if usingRNA:
            variantDict["RNASupport"][variantHash] = variantDict[variantHash]["RNASupport"]
    return variantDict
        
def collectTandemSNVSites(variantDict, maxDifferenceInPercentage = 0.10, maxFusionLength = 0):
    import operator
    variantHashList = list(variantDict["combined"].keys())
    variantHashList.sort(key = operator.itemgetter(1))
    variantHashList.sort(key = operator.itemgetter(0))
    tandemSiteTable = {}
    currentIndex = 0
    while currentIndex < len(variantHashList):
        currentHash = variantHashList[currentIndex]
        if currentHash[1] == 633353:
            print("something")
        currentContig, currentPosition, currentRef, currentAlt = currentHash
        currentTumorSupportingPercent = variantDict["combined"][currentHash].tumorSupportingPercent
        if len(currentRef) != 1 or len(currentAlt) != 1:  #catch for indels
            currentIndex += 1
            continue
        tandem = True
        nextIndex = currentIndex + 1
        tandemSites = []
        while tandem:
            if maxFusionLength and len(tandemSites) >= maxFusionLength + 1:
                currentIndex = nextIndex
                continue
            if nextIndex >= len(variantHashList):
                currentIndex = nextIndex
                break
            nextHash = variantHashList[nextIndex]
            nextContig, nextPosition, nextRef, nextAlt = nextHash
            nextTumorSupportingPercent = variantDict["combined"][nextHash].tumorSupportingPercent
            if len(nextRef) != 1 or len(nextAlt) != 1:
                nextIndex += 1
                continue
            if (currentContig == nextContig and nextPosition == currentPosition + 1) and (abs(currentTumorSupportingPercent - nextTumorSupportingPercent) <= maxDifferenceInPercentage):
                tandemSites.append(nextHash)
                tandem = True #not strictly necessary, but keeping to make it clear
                currentIndex = nextIndex
                currentTumorSupportingPercent = nextTumorSupportingPercent
                nextIndex += 1
                continue
            else:
                currentIndex = nextIndex
                tandem = False
                continue
        if tandemSites:
            tandemSiteTable[currentHash] = tandemSites
    return tandemSiteTable
            
def fuseTandemSites(tandemSiteHashDict, unfusedVariantDict):
    import variantDataHandler
    import scipy
    unfusedVariantDict["fused"] = {}
    unfusedVariantDict["fused"]["variants"] = unfusedVariantDict['combined'].copy()
    if "RNASupport" in unfusedVariantDict:
        usingRNA = True
    else:
        usingRNA = False
    if usingRNA:
        unfusedVariantDict["fused"]["RNASupport"] = unfusedVariantDict["RNASupport"].copy()
    sitesToFuse = list(tandemSiteHashDict.keys())
    for siteToFuse in sitesToFuse:
        newSiteData = unfusedVariantDict["combined"][siteToFuse] #using the existing site as a base
        if usingRNA:
            newSiteRNAData = unfusedVariantDict["RNASupport"][siteToFuse]
        tandemSites = tandemSiteHashDict[siteToFuse]
        for tandemSite in tandemSites:
            appendedSiteData = unfusedVariantDict['combined'][tandemSite]
            newSiteData.ref += appendedSiteData.ref
            newSiteData.alt += appendedSiteData.alt
            newSiteData.tumorDepth += appendedSiteData.tumorDepth  #we will average this out when all values have been added
            newSiteData.normalDepth += appendedSiteData.normalDepth
            newSiteData.tumorSupporting += appendedSiteData.tumorSupporting
            newSiteData.normalSupporting += appendedSiteData.normalSupporting
            if usingRNA:
                appendedSiteRNAData = unfusedVariantDict['RNASupport'][tandemSite]
                newSiteRNAData.score += appendedSiteRNAData.score
                newSiteRNAData.supportingReads += appendedSiteRNAData.supportingReads
                newSiteRNAData.totalDepth += appendedSiteRNAData.totalDepth
        newSiteData.tumorDepth = round(newSiteData.tumorDepth / len(tandemSites))
        newSiteData.tumorSupporting = round(newSiteData.tumorSupporting / len(tandemSites))
        newSiteData.normalDepth = round(newSiteData.normalDepth / len(tandemSites))
        newSiteData.normalSupporting = round(newSiteData.normalSupporting / len(tandemSites))
        if usingRNA:
            newSiteRNAData.score = round(newSiteRNAData.score / len(tandemSites))
            newSiteRNAData.supportingReads = round(newSiteRNAData.supportingReads / len(tandemSites))
            newSiteRNAData.totalDepth = round(newSiteRNAData.totalDepth / len(tandemSites))
            newSiteRNAData.oddsRatio, newSiteRNAData.pvalue = scipy.stats.fisher_exact([newSiteRNAData.supportingReads, newSiteRNAData.totalDepth],[newSiteData.tumorSupporting, newSiteData.tumorDepth])
        newSiteData.fusedSNV = True
        unfusedVariantDict["fused"]["variants"][siteToFuse] = newSiteData
        if usingRNA:
            unfusedVariantDict["fused"]["RNASupport"][siteToFuse] = newSiteRNAData
        for tandemSite in tandemSites:
            del unfusedVariantDict["fused"]["variants"][tandemSite]
            if usingRNA:
                del unfusedVariantDict["fused"]["RNASupport"][tandemSite]
    return unfusedVariantDict  #not sure if this is strictly necessary, as I think the function will be working on a reference to the original
        
def createOutputTextTable(variantDict):
    fusedDict = variantDict["fused"]
    headerLine = "\t".join(["contig", "position", "ref_allele", "alt_allele", "tumorSupporting", "tumorDepth", "normalSupporting", "normalDepth"])
    if "RNASupport" in variantDict:
        usingRNA = True
        headerLine += "\t" + "rnaSupport"
    else:
        usingRNA = False
    outputTable = [headerLine]
    for variant in list(fusedDict["variants"].keys()):
        outputLine = (str(fusedDict["variants"][variant]))
        if usingRNA:
            outputLine += "\t" + (str(fusedDict["RNASupport"][variant]))
        outputTable.append(outputLine)
    return outputTable

def moveRawVariantData(variantDict):
    keys = list(variantDict.keys())
    variantDict["rawVariants"] = {}
    for key in keys:
        if type(key) == tuple:
            variantDict["rawVariants"][key] = variantDict[key]
            del variantDict[key]
    return variantDict

def main():
    import pickle
    args = CheckArgs()
    variantDictFile = open(args.variantPickleFile, 'rb')
    variantDict = pickle.load(variantDictFile)
    variantDictFile.close()
    variantDict = makeCombinedAndRNAEntries(variantDict)
    tandemSiteTable = collectTandemSNVSites(variantDict, args.maxDifferencePercent, args.maxFusionLength)
    variantDict = fuseTandemSites(tandemSiteTable, variantDict)
    variantDict = moveRawVariantData(variantDict)
    if args.outputVariantFile.endswith(".pkl"):
        outputFile = open(args.outputVariantFile, 'wb')
        pickle.dump(variantDict, outputFile)
        outputFile.close()
    elif args.outputVariantFile.endswith(".txt"):
        outputFile = open(args.outputVariantFile, 'w')
        outputTextTable = createOutputTextTable(variantDict)
        for line in outputTextTable:
            print(line, file = outputFile)
        outputFile.close()            
    quit()
    
if __name__ == '__main__':
    main()
    