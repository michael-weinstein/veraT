#!/usr/bin/env python3

class VariantData(object):
    
    def __init__(self, contig, position, ref, alt, readDepth, supportingReads):
        self.contig = contig
        self.position = int(position)
        self.ref = ref
        self.alt = alt
        self.isIndel = len(self.ref) > 1 or len(self.alt) > 1 or "*" in alt
        self.depth = int(readDepth)
        self.supporting = int(supportingReads)
        if self.depth:
            self.supportingPercent = self.supporting / self.depth
        else:
            self.supportingPercent = 0
        self.hashValue = (self.contig, self.position, self.ref, self.alt)
        self.fusedSNV = False
    
    def __hash__(self):
        return hash(self.hashValue)
    
    def __eq__(self, other):
        if not type(other) == VariantData:
            raise TypeError("Can only compare VariantData types like this.  Passed value was %s type and had value %s" %(type(other), other))
        return self.hashValue == other.hashValue
    
    def __ne__(self, other):
        return not self.__eq__(other)
            
class SomaticVariantData(object):
    
    def __init__(self, contig, position, ref, alt, normalReadDepth, normalSupportingReads, tumorReadDepth, tumorSupportingReads):
        self.contig = contig
        self.position = int(position)
        self.ref = ref
        self.alt = alt
        self.isIndel = len(self.ref) > 1 or len(self.alt) > 1 or "*" in alt
        self.normalDepth = int(normalReadDepth)
        self.normalSupporting = int(normalSupportingReads)
        if self.normalDepth:
            self.normalSupportingPercent = self.normalSupporting / self.normalDepth
        else:
            self.normalSupportingPercent = 0
        self.tumorDepth = int(tumorReadDepth)
        self.tumorSupporting = int(tumorSupportingReads)
        self.readCountString = "%s,%s,%s,%s" %(self.tumorSupporting, self.tumorDepth, self.normalSupporting, self.normalDepth)
        if self.tumorDepth:
            self.tumorSupportingPercent = self.tumorSupporting / self.tumorDepth
        else:
            self.tumorSupportingPercent = 0
        self.pvalue = None
        self.hashValue = (self.contig, self.position, self.ref, self.alt)
        self.fusedSNV = False

    def fisherTest(self):
        if self.pvalue:
            return self.pvalue
        else:
            import scipy.stats
            '''
                    supporting  depth
            normal  counts      counts
            tumor   counts      counts
            '''
            self.oddsRatio, self.pvalue = scipy.stats.fisher_exact([[self.normalSupporting, self.normalDepth], [self.tumorSupporting, self.tumorDepth]])
            return self.pvalue
    
    def __hash__(self):
        return hash(self.hashValue)
    
    def __eq__(self, other):
        if not type(other) == VariantData:
            raise TypeError("Can only compare VariantData types like this.  Passed value was %s type and had value %s" %(type(other), other))
        return self.hashValue == other.hashValue
    
    def __ne__(self, other):
        return not self.__eq__(other)
    
    def __str__(self):
        printItems = [self.contig, self.position, self.ref, self.alt, self.tumorSupporting, self.tumorDepth, self.normalSupporting, self.normalDepth]
        printItems = [str(item) for item in printItems]
        return "\t".join(printItems)
    
    def oncotatorInputLine(self, delimiter = "\t"):
        endPosition = self.position + len(self.ref) - 1
        lineItems = [self.contig, self.position, endPosition, self.ref, self.alt, self.tumorDepth, self.tumorSupporting, self.tumorSupportingPercent, self.normalDepth]
        lineItems = [str(item) for item in lineItems]
        return delimiter.join(lineItems)
    
class RNASupportData(object):
    
    def __init__(self, score, supportingReads, totalDepth, pvalue = None, oddsRatio = None):
        self.score = score
        self.supportingReads = supportingReads
        self.totalDepth = totalDepth
        self.pvalue = pvalue
        self.oddsRatio = oddsRatio
    
    def __str__(self):
        printItems = [self.score, self.supportingReads, self.totalDepth, self.pvalue]
        printItems = [str(item) for item in printItems]
        return ",".join(printItems)

class ProteinChange(object):
    
    def __init__(self, gene, enst, changeTupleList):
        self.gene = gene
        self.enst = enst
        self.changes = changeTupleList
        self.polypeptideChange = len(changeTupleList) > 1
        
class OncotatorData(object):
    
    def __init__(self, gene, transcript, description, variantClassification, variantType, genomeChange, exon, cDNAChange, proteinChange):
        self.gene = gene
        self.transcript = transcript
        self.description = description
        self.variantClassification = variantClassification
        self.variantType = variantType
        self.genomeChange = genomeChange
        self.exon = exon
        self.cDNAChange = cDNAChange
        self.proteinChange = proteinChange
        
class NetMHCPrediction(object):
    
    def __init__(self, predictionLine):
        self.bindingFlag = None
        originalPredictionLine = predictionLine
        if type(predictionLine) == str:
            predictionLine = predictionLine.split()
            originalPredictionLine = predictionLine
        elif type(predictionLine) in (list, tuple):
            originalPredictionLine = predictionLine.copy()
        else:
            raise RuntimeError("Raw prediction line must be either string or list/tuple.  Got %s\n%s" %(type(predictionLine), predictionLine))
        if len(predictionLine) == 16:
            self.bindingFlag = predictionLine[-1]
            predictionLine = predictionLine[:14]
        elif len(predictionLine) == 14:
            pass
        else:
            raise RuntimeError("Got an invalid number of fields for prediction line %s" %originalPredictionLine)
        self.pos = predictionLine.pop(0)
        self.hla = predictionLine.pop(0).replace("HLA-","")
        self.peptide = predictionLine.pop(0)
        self.core = predictionLine.pop(0)
        self.offset = predictionLine.pop(0)
        self.ipos = int(predictionLine.pop(0))
        self.ilen = int(predictionLine.pop(0))
        self.dpos = int(predictionLine.pop(0))
        self.dlen = int(predictionLine.pop(0))
        self.icore = predictionLine.pop(0)
        self.identity = predictionLine.pop(0)
        self.log = float(predictionLine.pop(0))
        self.affinity = float(predictionLine.pop(0))
        self.rank = float(predictionLine.pop(0))
        
def sortVariantDataTuples(variantDataTupleList):
    if not type(variantDataTupleList) == list:
        raise RuntimeError("Variant data tuple list must be passed as a list type.")
    import operator
    def contigSortValue(rawContig, highestNumber = 99):
        import re
        digits = len(str(highestNumber))
        contig = re.sub("chr", "", rawContig.strip(), flags=re.IGNORECASE)
        sortOrder = "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y,M,MT".split(",")
        sortingTable = {}
        for index, contigName in enumerate(sortOrder):
            sortingTable[contigName] = str(index)
        if not contig in sortingTable:
            return rawContig
        else:
            returnValue = sortingTable[contig]
            if returnValue.isdigit():
                return returnValue.zfill(digits)
            else:
                return returnValue
    def contigSortValueFromSomaticTuple(somaticTuple, highestNumber = 99):
        return contigSortValue(somaticTuple[0], highestNumber)
    #goal is to have variants sorted first by contig in the standard order, then by position within the contig.  This will be done by first sorting on position and then on contig using a function that specifies the order.
    variantDataTupleList.sort(key = operator.itemgetter(1))  #sort first on position
    variantDataTupleList.sort(key = contigSortValueFromSomaticTuple)
