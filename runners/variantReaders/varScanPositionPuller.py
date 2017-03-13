#!/usr/bin/env python3

class CheckArgs():  #class that checks arguments and ultimately returns a validated set of arguments to the main program
    
    def __init__(self):
        import argparse
        import os
        parser = argparse.ArgumentParser()
        parser.add_argument("-f", "--file", help = "Variant call file path", required = True)
        parser.add_argument("-o", "--output", help = "Output pickle file name")
        rawArgs = parser.parse_args()
        file = rawArgs.file
        if not os.path.isfile(file):
            raise FileNotFoundError("Unable to find file %s" %file)
        self.file = file
        output = rawArgs.output
        if not output:
            output = file + ".positions"
        self.output = output
        
class VarScanColumnHeader(object):
    
    def __init__(self, rawLine):
        self.rawLine = rawLine.strip().strip("#")
        self.lineArray = self.rawLine.split("\t")
        self.fieldDict = {}
        self.fieldList = []
        specialFields = {"CHROM" : ["CHROMOSOME", "CONTIG", "CHR"],
                         "POSITION" : ["POS"],
                         "REF" : ["REFERENCE","REFERENCEALLELE","REFALLELE"],
                         "VAR" : ["ALT", "ALTERNATE", "ALTERNATIVE", "ALTALLELE", "ALTERNATIVEALLELE", "VARIANT", "VARIANTALLELE"]}
        for index, field in enumerate(self.lineArray):
            field = field.upper()
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
        self.header = header
        self.rawLine = rawLine.strip()
        self.lineArray = self.rawLine.split("\t")
        self.contig = self.lineArray[header.chrom]
        self.position = int(self.lineArray[header.pos])
        self.referenceAllele = self.lineArray[header.ref]
        self.end = self.position + len(self.referenceAllele)

if __name__ == '__main__':
    args = CheckArgs()
    varScanData = open(args.file, 'r')
    output = open(args.output, 'w')
    header = False
    for line in varScanData:
        if not header:
            header = VarScanColumnHeader(line)
        else:
            lineData = VarScanDataLine(line, header)
            print("%s\t%s\t%s" %(lineData.contig, lineData.position, lineData.end), file = output)
    varScanData.close()
    output.close()