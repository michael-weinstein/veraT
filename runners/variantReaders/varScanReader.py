#!/usr/bin/env python3

class CheckArgs():  #class that checks arguments and ultimately returns a validated set of arguments to the main program
    
    def __init__(self):
        import argparse
        import os
        parser = argparse.ArgumentParser()
        parser.add_argument("-f", "--file", help = "Variant call file path", required = True)
        parser.add_argument("-o", "--output", help = "Output pickle file name", required = True)
        parser.add_argument("-d", "--minDepth", help = "Minimum depth of coverage in both tumor and normal.", default = 10, type = int)
        rawArgs = parser.parse_args()
        file = rawArgs.file
        if not os.path.isfile(file):
            raise FileNotFoundError("Unable to find file %s" %file)
        self.file = file
        output = rawArgs.output
        if not output:
            output = self.file + "accepted.pkl"
        self.output = output
        self.minDepth = rawArgs.minDepth
        
class VarScanColumnHeader(object):  #need to use a preset column header here because the false positive filter
    
    def __init__(self):
        self.rawLine = "chrom,position,ref,var,normal_reads1,normal_reads2,normal_var_freq,normal_gt,tumor_reads1,tumor_reads2,tumor_var_freq,tumor_gt,somatic_status,variant_p_value,somatic_p_value,tumor_reads1_plus,tumor_reads1_minus,tumor_reads2_plus,tumor_reads2_minus,normal_reads1_plus,normal_reads1_minus,normal_reads2_plus,normal_reads2_minus"
        self.lineArray = self.rawLine.split(",")
        self.lineArray = [item.upper() for item in self.lineArray]
        self.fieldDict = {}
        self.fieldList = []
        specialFields = {"CHROM" : ["CONTIG", "CHROMOSOME", "CHR"],
                         "POSITION" : ["POS"],
                         "REF" : ["REF_ALLELE", "REFERENCE","REFERENCEALLELE","REFALLELE"],
                         "VAR" : ["ALT", "ALT_ALLELE", "ALTERNATE", "ALTERNATIVE", "ALTALLELE", "ALTERNATIVEALLELE", "VARIANT", "VARIANTALLELE"],
                         "TUMOR_READS1" : ["TUMORREFREADS", "TUMORREFCOUNT"],
                         "TUMOR_READS2" : ["TUMORALTREADS", "TUMORALTCOUNT"],
                         "NORMAL_READS1" : ["NORMALREFREADS", "NORMALREFCOUNT"],
                         "NORMAL_READS2" : ["NORMALALTREADS", "NORMALALTCOUNT"]}
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
                    print(self.fieldDict)
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
    
class VarScanDataLine(object):
    
    def __init__(self, rawLine, header):
        import variantDataHandler
        lineArray = rawLine.strip().split("\t")
        self.contig = lineArray[header.contig]
        self.position = int(lineArray[header.position])
        self.ref = lineArray[header.ref]
        self.alt = lineArray[header.alt]
        self.tumorRefCount = int(lineArray[header.tumorRefReads])
        self.tumorAltCount = int(lineArray[header.tumorAltReads])
        self.normalRefCount = int(lineArray[header.normalRefReads])
        self.normalAltCount = int(lineArray[header.normalAltReads])
        self.normalDepth = self.normalAltCount + self.normalRefCount
        self.tumorDepth = self.tumorAltCount + self.tumorRefCount
        self.somaticChange = variantDataHandler.SomaticVariantData(self.contig, self.position, self.ref, self.alt, self.normalDepth, self.normalAltCount, self.tumorDepth, self.tumorAltCount)

    
def analyzeVarScanSomatics(fileName, minDepth = 10, maxPValue = False):
    varScanData = open(fileName, 'r')
    somaticCollector = {}
    header = VarScanColumnHeader()
    for line in varScanData:
        line = line.strip()
        if not line:
            continue
        elif line.startswith("#"):
            continue
        else:
            data = VarScanDataLine(line, header)
            somaticChange = data.somaticChange
            if not somaticChange:
                continue
            if not (somaticChange.tumorDepth >= minDepth and somaticChange.normalDepth >= minDepth):
                continue
            somaticCollector[somaticChange.hashValue] = somaticChange
    varScanData.close()
    return somaticCollector

if __name__ == '__main__':
    args = CheckArgs()
    somaticList = analyzeVarScanSomatics(args.file, args.minDepth)
    print("Found %s somatic candidates" %len(somaticList))
    output = open(args.output, 'wb')
    import pickle
    pickle.dump(somaticList, output)
    output.close()