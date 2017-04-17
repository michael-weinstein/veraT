#!/usr/bin/env python3

class CheckArgs():  #class that checks arguments and ultimately returns a validated set of arguments to the main program
    
    def __init__(self):
        import argparse
        import os
        parser = argparse.ArgumentParser()
        parser.add_argument("-f", "--file", help = "Oncotator variant output file", required = True)
        parser.add_argument("-v", "--variantPickle", help = "Pickle containing variant data for sample", required = True)
        parser.add_argument("-o", "--output", help = "Output pickle file name")
        rawArgs = parser.parse_args()
        file = rawArgs.file
        if not os.path.isfile(file):
            raise FileNotFoundError("Unable to find oncotator variant file %s" %file)
        self.file = file
        variantPickle = rawArgs.variantPickle
        if not os.path.isfile(variantPickle):
            raise FileNotFoundError("Unable to find variant pickle file %s" %variantPickle)
        self.variantPickle = variantPickle
        output = rawArgs.output
        if not output:
            output = ".".join(self.file.split(".")[:-1]) + ".proteinChanges.pkl"
        self.output = output
        
class OncotatorColumnHeader(object):
    
    def __init__(self, rawLine):
        self.rawLine = rawLine.strip().strip("#")
        self.lineArray = self.rawLine.split("\t")
        self.fieldDict = {}
        self.fieldList = []
        specialFields = {"CHROMOSOME" : ["CONTIG", "CHR", "CHROM"],
                         "START_POSITION" : ["POSITION", "POS"],
                         "REFERENCE_ALLELE" : ["REF", "REFERENCE", "REFERENCEALLELE", "REFALLELE"],
                         "TUMOR_SEQ_ALLELE2" : ["ALT","ALTERNATE", "ALTERNATIVE", "ALTALLELE", "ALTERNATIVEALLELE", "VARIANT", "VARIANTALLELE"],  #if there are two variants for the position, we will turn this into a list
                         "VARIANT_CLASS" : ["VARIANTCLASS", "CLASS"],
                         "VARIANT_TYPE" : ["TYPE", "VARIANTTYPE"],
                         "PROTEIN_CHANGE": ["PROTEINCHANGE", "AMINOACIDCHANGE", "PEPTIDECHANGE"],
                         "ANNOTATION_TRANSCRIPT": ["ENST"],
                         "HUGO_SYMBOL" : ["GENE"]}
        for index, field in enumerate(self.lineArray):
            self.fieldDict[field] = index
            self.fieldList.append(field)
            if field.upper() in specialFields:
                for altName in specialFields[field.upper()]:
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
                    raise KeyError("%s is not a valid field in this dataset." %key)
                
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
    
class ProteinChangeHandler(object):
    
    def __init__(self, rawChange):
        rawChange = rawChange.replace("p.","")
        rawChange = rawChange.strip()
        if "*" in rawChange:
            self.nonsense = True
        else:
            self.nonsense = False
        if not "_" in rawChange:
            self.dipeptideChange = False
            self.reference = rawChange[0]
            self.variant = rawChange[-1]
            self.position = int(rawChange[1:-1])
            self.changes = [(self.reference, self.position, self.variant)]
            self.mutated = self.reference != self.variant
        else:
            import re
            self.dipeptideChange = True
            rawPositions = re.search("^\d+\_\d+", rawChange).group(0)
            rawPeptides = re.search("\D+\>\D+$", rawChange).group(0)
            self.positions = [int(number) for number in rawPositions.split("_")]
            self.reference, self.variant = rawPeptides.split(">")
            self.changes = []
            if self.reference == self.variant:
                self.mutated = False
            else:
                self.mutated = True
                startPosition = self.position[0]
                for i in range(0, len(self.reference)):
                    currentRef = self.reference[i]
                    currentPos = startPosition + i
                    currentAlt = self.variant[i]
                    if currentRef == currentAlt:
                        continue
                    self.changes.append((currentRef, currentPos, currentAlt))
                    
    def __bool__(self):
        return self.mutated
        
class OncotatorDataLine(object):
    
    def __init__(self, rawLine, header):
        self.header = header
        self.rawLine = rawLine.strip()
        self.lineArray = self.rawLine.split("\t")
        self.contig = self.lineArray[header.chrom]
        self.position = int(self.lineArray[header.pos])
        self.referenceAllele = self.lineArray[header.ref]
        self.altAllele = self.lineArray[header.alt]
        self.hashValue = (self.contig, self.position, self.referenceAllele, self.altAllele)
        self.proteinChange = self.lineArray[header.proteinchange]
        self.proteinChange = self.proteinChange.strip()
        if self.proteinChange:
            self.proteinChange = ProteinChangeHandler(self.proteinChange)
        self.enst = self.lineArray[header.enst]
        self.gene = self.lineArray[header.gene]
        
def analyzeOncotatorOutput(fileName, skipNonsense = True):
    import variantDataHandler
    oncotatorData = open(fileName, 'r')
    peptideChangeCollector = {}
    header = False
    for line in oncotatorData:
        line = line.strip()
        if not line:
            continue
        elif line.startswith("#"):
            continue
        elif not header:
            header = OncotatorColumnHeader(line)
            continue
        else:
            data = OncotatorDataLine(line, header)
            if not data.proteinChange:
                continue
            if skipNonsense and data.proteinChange.nonsense:
                continue
            peptideChangeCollector[data.hashValue] = variantDataHandler.ProteinChange(data.gene, data.enst, data.proteinChange.changes)
    oncotatorData.close()
    return peptideChangeCollector

if __name__ == '__main__':
    import pickle
    args = CheckArgs()
    peptideChanges = analyzeOncotatorOutput(args.file)
    print("Found %s peptide changes.  Adding to variant pickle." %len(peptideChanges))
    variantPickleFile = open(args.variantPickle, 'rb')
    variantPickle = pickle.load(variantPickleFile)
    variantPickleFile.close()
    variantPickle["fused"]["proteinChanges"] = peptideChanges
    output = open(args.output, 'wb')
    pickle.dump(variantPickle, output)
    output.close()