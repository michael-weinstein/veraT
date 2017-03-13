#!/usr/bin/env python3

class CheckArgs():  #class that checks arguments and ultimately returns a validated set of arguments to the main program
    
    def __init__(self):
        import argparse
        import os
        parser = argparse.ArgumentParser()
        parser.add_argument("-f", "--file", help = "Variant call file path", required = True)
        parser.add_argument("-t", "--tumor", help = "Tumor name", required = True)
        parser.add_argument("-n", "--normal", help = "Normal name", required = True)
        parser.add_argument("-o", "--output", help = "Output pickle file name")
        parser.add_argument("-p", "--maxPValue", help = "Maximum p-value in fisher test for variant.", default = 0.05, type = float)
        parser.add_argument("-d", "--minDepth", help = "Minimum depth of coverage in both tumor and normal.", default = 10, type = int)
        rawArgs = parser.parse_args()
        file = rawArgs.file
        if not os.path.isfile(file):
            raise FileNotFoundError("Unable to find file %s" %file)
        self.file = file
        self.tumor = rawArgs.tumor
        self.normal = rawArgs.normal
        output = rawArgs.output
        if not output:
            output = self.file + ".accepted.pkl"
        self.output = output
        self.maxPValue = rawArgs.maxPValue
        self.minDepth = rawArgs.minDepth
        
class VCFColumnHeader(object):
    
    def __init__(self, rawLine):
        self.rawLine = rawLine.strip().strip("#")
        self.lineArray = self.rawLine.split("\t")
        self.fieldDict = {}
        self.fieldList = []
        specialFields = {"CHROM" : ["CHROMOSOME", "CONTIG", "CHR"],
                         "POS" : ["POSITION"],
                         "REF" : ["REFERENCE","REFERENCEALLELE","REFALLELE"],
                         "ALT" : ["ALTERNATE", "ALTERNATIVE", "ALTALLELE", "ALTERNATIVEALLELE", "VARIANT", "VARIANTALLELE"],
                         "QUAL" : ["QUALITY"]}
        for index, field in enumerate(self.lineArray):
            self.fieldDict[field] = index
            self.fieldList.append(field)
            if field in specialFields:
                for altName in specialFields[field]:
                    self.fieldDict[altName] = index
        self.caseSensitives = self.getCaseCollisions()
        for index, field in enumerate(self.lineArray):
            if not field.upper() in self.caseSensitives:
                self.fieldDict[field.upper()] = index
        
    def getItemHandler(self, key):
        try:
            return self.fieldDict[key]
        except KeyError:
            if key.upper() in self.caseSensitives:
                raise KeyError("%s is in the list of fields, but is case sensitive due to a collision with another field.  This case did not match with a valid entry." %key)
            else:
                try:
                    return self.fieldDict[key.upper()]
                except KeyError:
                    raise KeyError("%s is not a valid field in this VCF." %key)
                
    def __getitem__(self, key):
        return self.getItemHandler(key)
    
    def __getattr__(self, attr):
        if attr == "fieldDict":
            return self.fieldDict
        elif attr == "fieldList":
            return self.fieldList
        else:
            return self.getItemHandler(attr)
        
    def getCaseCollisions(self):
        fieldList = self.fieldList
        upperFields = [item.upper() for item in fieldList]
        if len(upperFields) == len(set(upperFields)):
            return []
        caseSensitives = set()
        collector = set()
        for field in upperFields:
            if field in collector:
                caseSensitives.add(field)
            collector.add(field)
        return list(caseSensitives)
    
class VCFLineFormatData(object):
    
    def __init__(self, rawField):
        self.fieldList = rawField.strip().split(":")
        self.fieldList = [item.upper() for item in self.fieldList]
        specialFields = {"GT" : ["GENOTYPE"],
                         "AD" : ["ALLELEREADS", "ALLELICDEPTH", "ALLELEDEPTH"],
                         "DP" : ["READDEPTH","DEPTHOFCOVERAGE","COVERAGE"],
                         "GQ" : ["GENOTYPEQUALITY"],
                         "PL" : ["PHREDSCORES"]}
        self.fieldDict = {}
        for index, item in enumerate(self.fieldList):
            self.fieldDict[item] = index
            if item in specialFields:
                for altName in specialFields[item]:
                    self.fieldDict[altName] = index
    
    def getItemHandler(self, key):
        key = key.upper()
        return self.fieldDict[key]
    
    def __getitem__(self, key):
        return self.getItemHandler(key)
    
    def __getattr__(self, attr):
        if attr == "fieldList":
            return self.fieldList
        else:
            return self.getItemHandler(attr)
    
class VCFSampleData(object):
    
    def __init__(self, rawField, formatData, minimumDepth = 10):
        if rawField.strip().startswith("./."):
            self.called = False
        else:
            self.called = True
            self.analyzeData(rawField.strip(), formatData)
            
    def __bool__(self):
        return self.called
            
    def analyzeData(self, rawField, formatData):
        self.abnormal = False
        self.fieldSplit = rawField.strip().split(":")
        self.rawGenotype = self.fieldSplit[formatData.gt]
        self.rawAlleleDepth = self.fieldSplit[formatData.ad]
        genotypeProcessing = self.rawGenotype.split("/")
        genotypeProcessing = [int(item) for item in genotypeProcessing]
        genotypeProcessing = set(genotypeProcessing)
        self.genotype = sorted(list(genotypeProcessing))
        if len(self.genotype) == 2:
            self.heterozygous = True
            self.homozygous = False
        elif len(self.genotype) == 1:
            self.heterozygous = False
            self.homozygous = True
        else:
            self.heterozygous = False
            self.homozygous = False
            self.abnormal = True
        self.alleleDepths = self.rawAlleleDepth.split(",")
        self.alleleDepths = [int(item) for item in self.alleleDepths]
        self.depth = sum(self.alleleDepths)
        
        

class VCFDataLine(object):
    
    def __init__(self, rawLine, header):
        self.header = header
        self.rawLine = rawLine.strip()
        self.lineArray = self.rawLine.split("\t")
        self.contig = self.lineArray[header.chrom]
        self.position = int(self.lineArray[header.pos])
        self.referenceAllele = self.lineArray[header.ref]
        self.altAlleleList = self.lineArray[header.alt].split(",")
        self.alleleList = [self.referenceAllele] + self.altAlleleList
        formatData = VCFLineFormatData(self.lineArray[header.format])
        for i in range(9, len(self.lineArray)):
            self.lineArray[i] = VCFSampleData(self.lineArray[i], formatData)
        
    def checkForSomatic(self, normal, tumor):
        tumorAlleleSet = set(tumor.genotype)
        normalAlleleSet = set(normal.genotype)
        somaticMutation = tumorAlleleSet.difference(normalAlleleSet)
        return sorted(list(somaticMutation))
    
    def analyzeSomaticChanges(self, normalName, tumorName, checkDepth = 10, allowSomaticChangeToReference = False):
        import variantDataHandler
        normal = self.lineArray[self.header[normalName]]
        tumor = self.lineArray[self.header[tumorName]]
        if not (tumor.called and normal.called):
            return False
        if checkDepth:
            if not (tumor.depth >= checkDepth and normal.depth >= checkDepth):
                return False
        somaticVariants = self.checkForSomatic(normal, tumor)
        if not somaticVariants:
            return False
        else:
            variantDataList = []
            for alleleNumber in somaticVariants:
                if not allowSomaticChangeToReference and alleleNumber == 0:
                    continue
                variantAllele = self.alleleList[alleleNumber]
                variantData = variantDataHandler.SomaticVariantData(self.contig, self.position, self.referenceAllele, variantAllele, normal.depth, normal.alleleDepths[alleleNumber], tumor.depth, tumor.alleleDepths[alleleNumber])
                variantDataList.append(variantData)
        return variantDataList
    
def analyzeVCFSomatics(fileName, tumorName, normalName, maxPValue = 0.05, minDepth = 10):
    vcf = open(fileName, 'r')
    somaticCollector = {}
    for line in vcf:
        line = line.strip()
        if not line:
            continue
        elif line.startswith("##"):
            continue
        elif line.startswith("#"):
            header = VCFColumnHeader(line)
            continue
        else:
            if not header:
                raise RuntimeError("Hit what appears to be a line of data before finding a valid column header line. Line:\n%s" %line)
            data = VCFDataLine(line, header)
            somaticChanges = data.analyzeSomaticChanges(normalName, tumorName, minDepth)
            if not somaticChanges:
                continue
            for somaticChange in somaticChanges:
                if maxPValue:
                    if somaticChange.fisherTest() > maxPValue:
                        continue
                somaticCollector[somaticChange.hashValue] = somaticChange
    vcf.close()
    return somaticCollector

if __name__ == '__main__':
    args = CheckArgs()
    somaticList = analyzeVCFSomatics(args.file, args.tumor, args.normal, args.maxPValue, args.minDepth) #"/Users/michaelweinstein/sourceCodeRepository/testData/JMS.vcf", "JMS_scalp_022213", "JMS_Normal")
    print("Found %s somatic candidates" %len(somaticList))
    output = open(args.output, 'wb')
    import pickle
    pickle.dump(somaticList, output)
    output.close()