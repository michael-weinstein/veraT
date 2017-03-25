#!/usr/bin/env python3

class CheckArgs():  #class that checks arguments and ultimately returns a validated set of arguments to the main program
    
    def __init__(self):
        import argparse
        import os
        parser = argparse.ArgumentParser()
        parser.add_argument("-v", "--variants", help = "Pickles of somatic variant dictionaries. Can be comma-separated into source, SNP/INDEL, path.  source and SNP/INDEL are optional, but a source is needed if SNP/INDEL is given.", dest = "variantFiles", action = "append", required = True)
        parser.add_argument("-m", "--minHits", help = "Minimum hits required across the variant info files to be included in output", type = int)
        parser.add_argument("-x", "--maxHits", help = "Maximum hits allowed for a file to be included in the output", type = int)
        parser.add_argument("-o", "--output", help = "Output file name", required = True)
        rawArgs = parser.parse_args()
        variantFiles = rawArgs.variantFiles
        self.variantFiles = [VariantPickleFile(fileData) for fileData in variantFiles]
        self.sourceCount = countSources(self.variantFiles)
        minHits = rawArgs.minHits
        maxHits = rawArgs.maxHits
        if (type(minHits) == type(None) and type(maxHits) == type(None)):
            minHits = self.sourceCount - 1
        if type(minHits) == int:
            if minHits > self.sourceCount:
                raise ValueError("Cannot use a higher minimum hit requirement (%s) than there are sources (%s)." %(minHits, self.sourceCount))
            if minHits < 0:
                raise ValueError("Minimum hits must be zero or greater.  %s was given as the value." %minHits)
        if type(maxHits) == int:
            if maxHits < 1:
                raise ValueError("Maximum allowed hits must be a positive integer. %s was given as the value." %maxHits)
            if type(minHits) == int:
                if maxHits < minHits:
                    raise ValueError("Maximum allowed hits cannot be less than minimum required hits (this is impossible).  Given values: max=%s  min=%s" %(maxHits, minHits))
            if maxHits > sourceCount:
                raise ValueError("Cannot use a higher maximum hit requirement (%s) than there are sources (%s)." %(maxHits, self.sourceCount))
        self.minHits = minHits
        self.maxHits = maxHits
        output = rawArgs.output
        self.output = output
        
class VariantPickleFile(object):
    
    def __init__(self, rawFileInfo):
        import os
        self.rawFileInfo = rawFileInfo
        self.fileDataList = rawFileInfo.split(",")
        if len(self.fileDataList) == 1:
            self.filePath = self.fileDataList[0]
            self.source = None
            self.varType = None
        elif len(self.fileDataList) == 2:
            self.source, self.filePath = self.fileDataList
            self.varType = None
        elif len(self.fileDataList) == 3:
            self.source, self.varType, self.filePath = self.fileDataList
        else:
            raise RuntimeError("Invalid file info (too many fields) given in argument --variants %s" %rawFileInfo)
        if self.source:
            self.source = self.source.upper()
        if self.varType:
            self.varType = self.varType.upper()
            if not self.varType in ("SNP", "SNV", "INDEL"):
                raise RuntimeError("Invalid file type specified (must be SNP or INDEL) in --variant %s" %rawFileInfo)
            if self.varType == "SNV":
                self.varType = "SNP"                
        if not os.path.isfile(self.filePath):
            raise FileNotFoundError("Unable to find file %s specified in --variant argument %s" %(self.filePath, rawFileInfo))

def countSources(variantFiles):
    sourceSet = set()
    for file in variantFiles:
        sourceSet.add(file.source)
    count = len(sourceSet)
    if None in sourceSet:
        noneCount = 0
        for file in variantFiles:
            if type(file.source) == type(None):
                noneCount += 1
        count += noneCount - 1 #subtracting the one for the one file with a None source that contributed to the set
    return count

def dictionaryHitCount(dictionaryList):
    if type(dictionaryList) in (list, tuple):
        keyLists = [list(dictionary.keys()) for dictionary in dictionaryList]
    if type(dictionaryList) == dict:
        keyLists = []
        for key in list(dictionaryList.keys()):
            keyLists.append(list(dictionaryList[key].keys()))
    collector = {}
    for keyList in keyLists:
        for key in keyList:
            if not key in collector:
                collector[key] = 1
            else:
                collector[key] += 1
    return collector

def dictionaryHitFilter(dictionaryList, minimumHits = None, maximumHits = None):
    hitCounts = dictionaryHitCount(dictionaryList)
    accepted = []
    for key in list(hitCounts.keys()):
        count = hitCounts[key]
        if not type(minimumHits) == type(None):
            if count < minimumHits:
                continue
        if not type(maximumHits) == type(None):
            if count > maximumHits:
                continue
        accepted.append(key)
    return accepted

def sortVariantDataTuples(variantDataTupleList):
    import operator
    import variantSupport
    #goal is to have variants sorted first by contig in the standard order, then by position within the contig.  This will be done by first sorting on position and then on contig using a function that specifies the order.
    variantDataTupleList.sort(key = operator.itemgetter(1))
    variantDataTupleList.sort(key = variantSupport.contigSortValueFromSomaticTuple)
    
def combineIndelAndSNP(indelDict, snpDict, priority = "SNP"):
    if priority.upper() == "SNP":
        newDict = indelDict.copy()
        newDict.update(snpDict)
    elif priority.upper() == "INDEL":
        newDict = snpDict.copy()
        newDict.update(indelDict)
    else:
        raise RuntimeError("Invalid priority set; must be either SNP or INDEL.  %s was given" %priority)
    return newDict
    
def sortVariantFiles(variantFiles):
    import os
    variantFileDict = {}
    counter = 0
    for file in variantFiles:
        source = file.source
        if type(source) == type(None):
            source = file.filePath.split(os.sep)[-1]
            counter += 1
        varType = file.varType
        path = file.filePath
        if not source in variantFileDict:
            variantFileDict[source] = {}
        if varType in variantFileDict[source]:
            raise RuntimeError("File collision involving %s\nDict: %s" %file.rawFileInfo)
        variantFileDict[source][varType] = path
    for source in list(variantFileDict.keys()):
        if len(variantFileDict[source]) == 1:
            assert set(variantFileDict[source]) == {None}, "Got a source group with unmated INDEL or SNP for %s\nDict: %s" %(source, variantFileDict[source].keys())
        elif len(variantFileDict[source]) == 2:
            assert set(variantFileDict[source].keys()) == {"INDEL", "SNP"}, "Got a source group with unmated INDEL or SNP for %s\nDict: %s" %(source, variantFileDict[source])
        else:
            raise RuntimeError("Got a source with inappropriate number of variant types %s\nDict: %s" %(source, variantFileDict[source]))
    return variantFileDict

def createOutputTextTable(sortedAcceptedVariantInfoTuples, variantDicts):
    sources = sorted(list(variantDicts.keys()))
    headerLine = ["contig", "position", "ref_allele", "alt_allele", "CountMethod"] + sources
    outputTable = [headerLine]
    for variant in sortedAcceptedVariantInfoTuples:
        contig, position, ref, alt = variant
        outputLine = [contig, position, ref, alt]
        countData = []
        hits = 0
        for source in sources:
            if variant in variantDicts[source]:
                countData.append(variantDicts[source][variant].readCountString)
                hits += 1
            else:
                countData.append("NA")
        outputLine = outputLine + [hits] + countData
        outputLine = (str(field) for field in outputLine)
        outputTable.append(outputLine)
    return outputTable

def roundedAverage(data):
    return round(sum(data)/len(data))

def createOutputPickleTable(sortedAcceptedVariantInfoTuples, variantDicts):
    import variantDataHandler
    sources = list(variantDicts.keys())
    genericEntryDict = {"hits" : 0}
    for source in sources:
        genericEntryDict[source] = False
    outputTable = {}
    for variant in sortedAcceptedVariantInfoTuples:
        normalDepth = []
        normalSupporting = []
        tumorDepth = []
        tumorSupporting = []
        outputTable[variant] = genericEntryDict.copy()
        for source in sources:
            if variant in variantDicts[source]:
                outputTable[variant][source] = variantDicts[source][variant]
                outputTable[variant]["hits"] += 1
                normalDepth.append(variantDicts[source][variant].normalDepth)
                normalSupporting.append(variantDicts[source][variant].normalSupporting)
                tumorDepth.append(variantDicts[source][variant].tumorDepth)
                tumorSupporting.append(variantDicts[source][variant].tumorSupporting)
        normalDepth = roundedAverage(normalDepth)
        normalSupporting = roundedAverage(normalSupporting)
        tumorDepth = roundedAverage(tumorDepth)
        tumorSupporting = roundedAverage(tumorSupporting)
        contig, position, ref, alt = variant
        outputTable[variant]["combined"] = variantDataHandler.SomaticVariantData(contig, position, ref, alt, normalDepth, normalSupporting, tumorDepth, tumorSupporting)
        pvalue = outputTable[variant]["combined"].fisherTest()
    return outputTable

def main():
    import pickle
    args = CheckArgs()
    variantFileDict = sortVariantFiles(args.variantFiles)
    variantDicts = {}
    for source in list(variantFileDict.keys()):
        if None in variantFileDict[source]:
            file = open(variantFileDict[source][None], 'rb')
            variantDicts[source] = pickle.load(file)
            file.close()
        elif set(variantFileDict[source].keys()) == {"SNP", "INDEL"}:
            snpFile = open(variantFileDict[source]["SNP"], 'rb')
            indelFile = open(variantFileDict[source]["INDEL"], 'rb')
            snpData = pickle.load(snpFile)
            indelData = pickle.load(indelFile)
            snpFile.close()
            indelFile.close()
            variantDicts[source] = combineIndelAndSNP(indelData, snpData)
        else:
            raise RuntimeError("Got a source with inappropriate variant types %s\nDict: %s" %(source, variantFileDict[source]))
    filteredVariants = dictionaryHitFilter(variantDicts, args.minHits, args.maxHits)
    sortVariantDataTuples(filteredVariants)
    if args.output.upper().endswith(".PKL"):
        outputTable = createOutputPickleTable(filteredVariants, variantDicts)
        outputFile = open(args.output, 'wb')
        pickle.dump(outputTable, outputFile)
        outputFile.close()
    else:
        outputTable = createOutputTextTable(filteredVariants, variantDicts)
        outputFile = open(args.output, 'w')
        for line in outputTable:
            print("\t".join(line), file = outputFile)
        outputFile.close()
    quit()
    
if __name__ == '__main__':
    main()
    