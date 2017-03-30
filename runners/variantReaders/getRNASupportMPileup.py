#!/usr/bin/env python3

class CheckArgs():  #class that checks arguments and ultimately returns a validated set of arguments to the main program
    
    def __init__(self):
        import argparse
        import os
        parser = argparse.ArgumentParser()
        parser.add_argument("-f", "--mpileupFile", help = "RNA mPileup", required = True)
        parser.add_argument("-s", "--somaticVariants", help = "Pickle containing somatic variant analysis from DNA", required = True)
        parser.add_argument("-o", "--output", help = "Output pickle file name")
        parser.add_argument("-m", "--minDiff", help = "Minimum percent difference in expression vs. DNA mutant/wild-type ratios to consider worth scoring", default = 10, type = int)
        rawArgs = parser.parse_args()
        mpileupFile = rawArgs.mpileupFile
        if not os.path.isfile(mpileupFile):
            raise FileNotFoundError("Unable to find file %s" %mpileupFile)
        self.mpileupFile = mpileupFile
        somaticVariants = rawArgs.somaticVariants
        if not os.path.isfile(somaticVariants):
            raise FileNotFoundError("Unable to find file %s" %somaticVariants)
        self.somaticVariants = somaticVariants
        output = rawArgs.output
        if not output:
            output = self.mpileupFile + ".rnaSupport.pkl"
        self.output = output
        self.minDiff = rawArgs.minDiff

def checkMPileupForRNASupport(mpileupFile, somaticVariants, somaticVariantTable,  minDiff):
    import variantDataHandler
    import scipy.stats
    import re
    nonBaseRegex = re.compile("[^ATGC]", re.IGNORECASE)  #compiling this regex now so it only needs to be compiled once
    plusBaseRegex = re.compile("\+\d+[ATGC]+", re.IGNORECASE)
    if mpileupFile.endswith(".gz"):
        import gzip
        mpileup = gzip.open(mpileupFile, 'rt')
    else:
        mpileup = open(mpileupFile, 'r')
    supportData = {}
    lociOfInterest = [(item[0], str(item[1])) for item in somaticVariants]  #using a string of the position to improve performance, will convert hundreds of times now instead of converting millions of times during mpileup reading
    for line in mpileup:
        line = line.strip()
        if not line:
            continue
        elif line.startswith("#"):
            continue
        else:
            try:
                contig, position, ref, depth, reads, quality = line.split()
            except ValueError:  #when lines with no reads are included, we get 4 columns (reads and qualities are missing)
                continue
            locus = (contig, position)
            if not locus in lociOfInterest:
                continue  #if we read past this, we know we are looking at a locus with a somatic change
            position = int(position)
            totalDepthRNA = int(depth)
            reads = reads.upper().replace(".", ref).replace(",", ref)
            reads = re.sub(plusBaseRegex, "", reads)  #get rid of any plus (digit)(base) reads here, otherwise the base will remain as an alt
            reads = re.sub(nonBaseRegex, "", reads)  #get rid of anything not a base (we have already replaced the periods and commas)
            foundHashesAtSite = []
            for variant in somaticVariants:
                if contig == variant[0] and position == variant[1]:
                    foundHashesAtSite.append(variant)
            for foundHash in foundHashesAtSite:
                if not (len(foundHash[2]) == 1 and len(foundHash[3]) == 1): #indel catcher
                    continue
                altAllele = foundHash[3]
                supportingDepthRNA = reads.count(altAllele)
                totalDepthDNA = somaticVariantTable[foundHash]["combined"].tumorDepth
                supportingDepthDNA = somaticVariantTable[foundHash]["combined"].tumorSupporting
                expressionDNARatio = (supportingDepthRNA/totalDepthRNA) / (supportingDepthDNA/totalDepthDNA)
                if not supportingDepthRNA:
                    supportData[foundHash] = variantDataHandler.RNASupportData(1, supportingDepthRNA, totalDepthRNA, None)
                    continue
                if totalDepthRNA < 10 or supportingDepthRNA < 3:
                    supportData[foundHash] = variantDataHandler.RNASupportData(2, supportingDepthRNA, totalDepthRNA, None)
                    continue
                if expressionDNARatio > (1 - minDiff/100) and expressionDNARatio < (1 + minDiff/100): #anything not greater or less than minDiff percent off expected will be called as not significantly different
                    supportData[foundHash] = variantDataHandler.RNASupportData(4, supportingDepthRNA, totalDepthRNA, None)
                    continue
                oddsRatio, pvalue = scipy.stats.fisher_exact([[supportingDepthRNA, totalDepthRNA], [supportingDepthDNA, totalDepthDNA]])
                if pvalue > 0.05:
                    supportData[foundHash] = variantDataHandler.RNASupportData(4, supportingDepthRNA, totalDepthRNA, pvalue)
                elif expressionDNARatio > 1:
                    supportData[foundHash] = variantDataHandler.RNASupportData(5, supportingDepthRNA, totalDepthRNA, pvalue)
                else:
                    supportData[foundHash] = variantDataHandler.RNASupportData(3, supportingDepthRNA, totalDepthRNA, pvalue)
    mpileup.close()
    return supportData
               
def createOutputTextTable(sortedAcceptedVariantInfoTuples, variantDicts):
    sources = sorted(list(variantDicts[sortedAcceptedVariantInfoTuples[0]].keys()))
    sources.remove("hits")
    sources.remove("combined")
    sources.remove("RNASupport")
    headerLine = ["contig", "position", "ref_allele", "alt_allele", "CountMethod"] + sources + ["RNASupport"]
    outputTable = [headerLine]
    for variant in sortedAcceptedVariantInfoTuples:
        contig, position, ref, alt = variant
        outputLine = [contig, position, ref, alt]
        countData = []
        hits = variantDicts[variant]["hits"]
        for source in sources:
            if source in variantDicts[variant] and variantDicts[variant][source]:
                countData.append(variantDicts[variant][source].readCountString)
            else:
                countData.append("NA")
        countData.append(str(variantDicts[variant]["RNASupport"]))
        outputLine = outputLine + [hits] + countData
        outputLine = (str(field) for field in outputLine)
        outputTable.append(outputLine)
    return outputTable

def sortVariantDataTuples(variantDataTupleList):
    import operator
    import variantSupport
    #goal is to have variants sorted first by contig in the standard order, then by position within the contig.  This will be done by first sorting on position and then on contig using a function that specifies the order.
    variantDataTupleList.sort(key = operator.itemgetter(1))
    variantDataTupleList.sort(key = variantSupport.contigSortValueFromSomaticTuple)
    
def main():
    args = CheckArgs()
    import pickle
    import variantDataHandler
    somaticsFile = open(args.somaticVariants, 'rb')
    somaticVariantTable = pickle.load(somaticsFile)
    somaticsFile.close()
    somaticVariants = list(somaticVariantTable.keys())
    for key in somaticVariants:
        somaticVariantTable[key]["RNASupport"] = variantDataHandler.RNASupportData(0, 0, 0)
    rnaSupportTable = checkMPileupForRNASupport(args.mpileupFile, somaticVariants, somaticVariantTable, args.minDiff)
    for key in list(rnaSupportTable.keys()):
        somaticVariantTable[key]["RNASupport"] = rnaSupportTable[key]
    sortVariantDataTuples(somaticVariants)
    if args.output.upper().endswith(".PKL"):
        outputFile = open(args.output, 'wb')
        pickle.dump(somaticVariantTable, outputFile)
        outputFile.close()
    else:
        outputTable = createOutputTextTable(somaticVariants, somaticVariantTable)
        outputFile = open(args.output, 'w')
        for line in outputTable:
            print("\t".join(line), file = outputFile)
        outputFile.close()
    quit()    
if __name__ == '__main__':
    main()
    