#!/usr/bin/env python3

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
    