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
    
class RNASupportData(object):
    
    def __init__(self, score, supportingReads, totalDepth, pvalue = None):
        self.score = score
        self.supportingReads = supportingReads
        self.totalDepth = totalDepth
        self.pvalue = pvalue
    
    def __str__(self):
        printItems = [self.score, self.supportingReads, self.totalDepth, self.pvalue]
        printItems = [str(item) for item in printItems]
        return "\t".join(printItems)
