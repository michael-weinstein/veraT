#!/usr/bin/env python3

class CheckArgs():  #class that checks arguments and ultimately returns a validated set of arguments to the main program
    
    def __init__(self):
        import argparse
        import os
        parser = argparse.ArgumentParser()
        parser.add_argument("-v", "--vcf", help = "RNA Variant call file path", required = True)
        parser.add_argument("-s", "--somaticVariants", help = "Pickle containing somatic variant analysis from DNA", required = True)
        parser.add_argument("-r", "--rna", help = "Tumor RNA variant column name", required = True)
        parser.add_argument("-o", "--output", help = "Output pickle file name")
        parser.add_argument("-m", "--minDiff", help = "Minimum percent difference in expression vs. DNA mutant/wild-type ratios to consider worth scoring", default = 10, type = int)
        rawArgs = parser.parse_args()
        vcf = rawArgs.vcf
        if not os.path.isfile(vcf):
            raise FileNotFoundError("Unable to find file %s" %vcf)
        self.vcf = vcf
        somaticVariants = rawArgs.somaticVariants
        if not os.path.isfile(somaticVariants):
            raise FileNotFoundError("Unable to find file %s" %somaticVariants)
        self.somaticVariants = somaticVariants
        self.rna = rawArgs.rna
        output = rawArgs.output
        if not output:
            output = self.vcf + ".rnaSupport.pkl"
        self.output = output
        self.minDiff = rawArgs.minDiff

def checkVCFForRNASupport(vcfPath, somaticVariants, somaticVariantTable, rnaSampleName, minDiff):
    import vcfReader
    import variantDataHandler
    import scipy.stats
    vcf = open(vcfPath, 'r')
    supportData = {}
    for line in vcf:
        line = line.strip()
        if not line:
            continue
        elif line.startswith("##"):
            continue
        elif line.startswith("#"):
            header = vcfReader.VCFColumnHeader(line)
            continue
        else:
            if not header:
                raise RuntimeError("Hit what appears to be a line of data before finding a valid column header line. Line:\n%s" %line)
            data = vcfReader.VCFDataLine(line, header, somaticVariants)
            if not data.inHashList:
                continue
            for foundHash in data.inHashList:
                if not data.lineArray[header[rnaSampleName]].called:
                    continue
                altAllele = foundHash[3]
                totalDepthRNA = data.lineArray[header[rnaSampleName]].depth
                supportingDepthRNA = data.lineArray[header[rnaSampleName]].alleleDepthTable[altAllele]
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
    vcf.close()
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
    rnaSupportTable = checkVCFForRNASupport(args.vcf, somaticVariants, somaticVariantTable, args.rna, args.minDiff)
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
    