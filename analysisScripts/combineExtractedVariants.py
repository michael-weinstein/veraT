#!/usr/bin/env python3

class CheckArgs(object):
   
   def __init__(self):
         import argparse
         import os
         import sys
         parser = argparse.ArgumentParser()
         parser.add_argument("-o", "--outputFile", help = "Output file.")
         parser.add_argument("-n", "--normal", help = "Require at least n reads supporting a variant for it to be emitted")
         parser.add_argument("-v", "--verbose", help = "Verbose mode", action = 'store_true')
         parser.add_argument("-t", "--tumor", help = "Filename for output of target list.  This puts the program into tumor sample mode.")
         rawArgs = parser.parse_args()
         if rawArgs.tumor:
            if os.path.isfile(rawArgs.tumor):
                self.tumor = rawArgs.tumor
            else:
                raise RuntimeError("Tumor file not found: " + rawArgs.tumor)
         else:
            raise RuntimeError("No tumor file specified.")
         if rawArgs.normal:
            if os.path.isfile(rawArgs.normal):
                self.normal = rawArgs.normal
            else:
                raise RuntimeError("Normal file not found: " + rawArgs.normal)
         else:
            raise RuntimeError("No normal file specified.")
         self.outputFile = rawArgs.outputFile
         if self.outputFile == self.tumor or self.outputFile == self.normal:
            raise RuntimeError("Error: Input and output files cannot be the same.")
         self.verbose = rawArgs.verbose
         
class VariantDataLine(object):
   
   def __init__(self, line, delimiter = "\t"):
      self.lineList = line.split("\t")
      self.mutantReadCount = int(self.lineList[0])
      self.coverage = int(self.lineList[1])
      self.mutantPercent = float(self.lineList[2])
      self.contig = self.lineList[3]
      self.position = int(self.lineList[4])
      self.locus = (self.contig, self.position)
      self.referenceBase = self.lineList[5]
      self.variant = self.lineList[6]
      self.variantType = self.lineList[7]
      self.strandDataString = self.lineList[8]
      workingStrand = self.strandDataString.replace("+","")
      workingStrand = workingStrand.split("-")
      self.plusReads, self.minusReads = [int(number) for number in workingStrand]
      del workingStrand
      self.rawPileupLine = self.lineList[9]
      self.rawPileupList = self.rawPileupLine.split("|")
      self.isIndel = False
      self.indelSeq = self.variant
      if self.variant[0] in "+-":
         self.isIndel = True
         self.analyzeIndel()
         
   def analyzeIndel(self):
      import re
      regex = re.compile("([+-])(\d+)(.+)")    #sorry
      result = re.match(regex, self.variant)
      self.indelType = result.group(1)
      self.indelLength = int(result.group(2))
      self.indelSeq = result.group(3)
      self.variant = ""
      if self.indelType == "+":
         self.variant += "I-"
      elif self.indelType == "-":
         self.variant += "D-"
      self.variant += str(self.indelLength)
      
class UncoveredDataLine(VariantDataLine):
   
   def __init__(self, locus, variant, variantType, delimiter = "\t"): #should take locus as an iterable of contig and position
      self.contig = locus[0]
      self.position = int(locus[1])
      self.mutantReadCount = 0
      self.coverage = 0
      self.mutantPercent = 0
      self.locus = (self.contig, self.position)
      self.referenceBase = "N"
      self.variant = "Unread"
      self.variantType = "U"
      self.strandDataString = "+0-0"
      self.plusReads = 0
      self.minusReads = 0
      self.rawPileupLine = self.contig + "|" + str(self.position) + "|N|0||"
      self.rawPileupList = self.rawPileupLine.split("|")
      self.lineList = [self.mutantReadCount, self.coverage, self.mutantPercent, self.contig, self.position, self.referenceBase, self.variant, self.variantType, self.strandDataString, self.contig + "|" + str(self.position) + "|N|0||"]

      
   def __nonzero__(self):  #an uncovered line class will always evaluate to false.  This should be convenient, but might cause issues depending on handling.
      return False
   
class UnsupportedDataLine(VariantDataLine):
   
   def __init__(self, line, variant, variantType, delimiter = "\t"):
      self.lineList = line.lineList
      self.mutantReadCount = 0
      self.coverage = int(self.lineList[1])
      self.mutantPercent = 0.0
      self.contig = self.lineList[3]
      self.position = int(self.lineList[4])
      self.locus = (self.contig, self.position)
      self.referenceBase = self.lineList[5]
      self.variant = variant  #replace with the one we are looking at which was not supported
      self.variantType = variantType
      self.strandDataString = "+0-0"
      workingStrand = self.strandDataString.replace("+","")
      workingStrand = workingStrand.split("-")
      self.plusReads, self.minusReads = [int(number) for number in workingStrand]
      del workingStrand
      self.rawPileupLine = self.lineList[9]
      self.rawPileupList = self.rawPileupLine.split("|")
      self.isIndel = False
      self.lineList = [self.mutantReadCount, self.coverage, self.mutantPercent, self.contig, self.position, self.referenceBase, self.variant, self.variantType, self.strandDataString, self.rawPileupLine]
      self.indelSeq = self.variant
      if self.variant[0] in "+-":
         self.isIndel = True
         self.analyzeIndel()
   
class ExtractedVariantFile(object):
   
   def __init__(self, fileHandle, dataType = "input"):
      self.fileHandle = fileHandle
      self.currentLine = False
      self.nextLine = False
      self.completed = False
      self.dataType = dataType
      self.buildContigTuple()
   
   def buildContigTuple(self):
      if args.verbose:
         print("Building contig list for the %s file." %(self.dataType))
      self.contigTuple = list()
      line = self.fileHandle.readline().strip()
      while line:
         if line.startswith("#"):
            line = self.fileHandle.readline().strip()
            continue
         contig = line.split("\t")[3]
         if not contig in self.contigTuple:
            self.contigTuple.append(contig)
         line = self.fileHandle.readline().strip()
      self.fileHandle.seek(0)
      self.contigTuple = tuple(self.contigTuple)
      
   def subsequentLineMatches(self):
      self.nextLine = False
      nextLine = self.fileHandle.readline().strip()
      if not nextLine:
         self.completed = True
         return False
      while nextLine.startswith("#"):
         nextLine = self.fileHandle.readline().strip()
         if not nextLine:
            self.completed = True
            return False
      self.nextLine = VariantDataLine(nextLine)
      if self.nextLine.locus == self.currentLine.locus:
         return True
      else:
         return False
      
   def loadNextLine(self):
      if self.nextLine:
         self.currentLine = self.nextLine
         self.nextLine = False
      else:
         currentLine = self.fileHandle.readline().strip()
         if not currentLine:
            self.completed = True
            return False
         while currentLine.startswith("#"):
            currentLine = self.fileHandle.readline().strip()
            if not currentLine:
               self.completed = True
               return False
         self.currentLine = VariantDataLine(currentLine)
      
   def searchForMatchedLocus(self, searchLocus, tumorContigList):
      if self.currentLine and self.currentLine.locus == searchLocus:  #woohoo, already there!
         return True
      searchContig, searchPosition = searchLocus
      while not self.currentLine and not self.completed:  #handler for when we don't have a currentLine loaded already
         self.loadNextLine()
      if self.completed:
         return False
      if self.currentLine.locus == searchLocus:
         return True
      searchContigIndex = tumorContigList.index(searchLocus[0])
      currentContigIndex = tumorContigList.index(self.currentLine.contig)
      while currentContigIndex < searchContigIndex and not self.completed:  #race ahead to the contig of interest if not already there
         self.loadNextLine()
         currentContigIndex = tumorContigList.index(self.currentLine.contig)
      if self.completed:
         return False
      if currentContigIndex > searchContigIndex:  #went past the contig of interest without finding it
         return False
      if self.currentLine.position > searchPosition:  #if we have already passed the point of interest
         return False
      while not self.currentLine.locus == searchLocus and not self.completed:  #keep running this until we hit the end of file or our match (traps for other stop conditions will be included in the loop)
         self.loadNextLine()  #we can do this because we already know we are not on the correct locus
         while not self.currentLine or self.completed:  #catch for completion or empty lines
            if self.completed:
               return False
            self.loadNextLine
         if self.currentLine.locus == searchLocus:  #woohoo, found it
            break
         if tumorContigList.index(self.currentLine.contig) > searchContigIndex: #went to a contig past the one of interest
            return False
         if self.currentLine.position > searchPosition:  #went past the position on the same contig without finding it
            return False
      if self.completed:  #in theory, this should have dropped out during the loop
         return False
      if self.currentLine.locus == searchLocus:
         return True
   
   def getMatchedLocusTable(self, locus, tumorContigList):
      if self.searchForMatchedLocus(locus, tumorContigList):
         locusVarTable = {self.currentLine.variant : self.currentLine}
         while self.subsequentLineMatches():
            self.loadNextLine()
            locusVarTable[self.currentLine.variant] = self.currentLine
         return locusVarTable
      else:
         return False
   
   def getNextLocusTable(self):
      if self.completed:
         return False
      self.loadNextLine()
      locusVarTable = {self.currentLine.variant : self.currentLine}
      while self.subsequentLineMatches():
         self.loadNextLine()
         locusVarTable[self.currentLine.variant] = self.currentLine
      return locusVarTable
      
class OutputLine(object):
   
   def __init__(self, tumorVariantLine, normalVariantLine, delimiter = "\t"):
      self.tumor = tumorVariantLine
      self.normal = normalVariantLine
      self.locus = tumorVariantLine.contig + ":" + str(tumorVariantLine.position)
      self.typeCodeTranslation = {"S" : "SUBSTITUTION",
                             "I" : "INDEL"}
      self.variantDataColumns = [self.locus, self.tumor.referenceBase, self.tumor.variant, self.typeCodeTranslation[self.tumor.variantType]]
      self.readDataColums = [self.tumor.mutantReadCount, self.tumor.coverage, self.normal.mutantReadCount, self.normal.coverage, self.tumor.strandDataString, self.normal.strandDataString, self.tumor.indelSeq]
      self.outputList = self.variantDataColumns + self.readDataColums + self.tumor.lineList + self.normal.lineList
      self.outputStringList = [str(item) for item in self.outputList]
      self.outputString = delimiter.join(self.outputStringList)
      
   def __str__(self):
      return self.outputString
   
def validContigLists(tumorContigTuple, normalContigTuple):   #checking to make sure that we don't have contigs out of order
   if tumorContigTuple == normalContigTuple:
      return True
   tumorContigSet = set(tumorContigTuple)
   normalContigSet = set(normalContigTuple)
   absentInTumor = normalContigSet - tumorContigSet
   absentInNormal = tumorContigSet - normalContigSet
   if not absentInTumor and not absentInNormal:  #this would mean that we have the same contigs in a different order
      return False  #DO NOT WANT
   removeFromTumorIndices = []
   removeFromNormalIndices = []
   tumorContigList = list(tumorContigTuple)
   normalContigList = list(normalContigTuple)
   for contig in absentInTumor:
      removeFromNormalIndices.append(normalContigList.index(contig))
   for contig in absentInNormal:
      removeFromTumorIndices.append(tumorContigList.index(contig))
   removeFromNormalIndices.sort(reverse = True)
   removeFromTumorIndices.sort(reverse = True)
   for index in removeFromNormalIndices:  #removals must go in reverse order (how we sorted the list) so we don't change the index of something before we delete it
      del normalContigList[index]
   for index in removeFromTumorIndices:
      del tumorContigList[index]
   if normalContigList == tumorContigList:
      return True
   else:
      return False

if __name__=="__main__":
   titleLine = "\t".join(["#locus", "ref", "var", "varType", "tumorVarSupport", "tumorCoverage", "normalVarSupport", "normalCoverage", "tumorStrands", "normalStrands", "varSeq", "tumorVarSup", "tumorCov", "tumorMutPercent", "tumorContig", "tumorPosition", "tumorRef", "tumorVar", "tumorVarType", "tumorStrand", "tumorPileup", "normalVarSup", "normalCov", "normalMutPercent", "normalContig", "normalPosition", "normalRef", "normalVar", "normalVarType", "normalStrand", "normalPileup"])
   import datetime
   startTime = datetime.datetime.now()
   global args
   args = CheckArgs()
   tumorFile = open(args.tumor, 'r')
   normalFile = open(args.normal, 'r')
   tumor = ExtractedVariantFile(tumorFile, "tumor")
   normal = ExtractedVariantFile(normalFile, "normal")
   outputFile = open(args.outputFile, 'w')
   print(titleLine, file = outputFile)
   if not validContigLists(tumor.contigTuple, normal.contigTuple):  #input validation is magic and lets us confirm that the contigs are not in some kind of different order between the two sample sets
      raise RuntimeError("Tumor and normal files do not appear to have been sorted in the same order.\nTumor contigs: %s\nNormal contigs: %s" %(tumor.contigTuple, normal.contigTuple))
   locusCount = 0
   while not tumor.completed:
      if args.verbose:
         if locusCount % 100 == 0:
            print("Comparing locus %s" %(locusCount), end = '\r')
      locusCount += 1
      tumorLocusVariantTable = tumor.getNextLocusTable()
      if not tumorLocusVariantTable:
         break
      matchedNormalVariantTable = normal.getMatchedLocusTable(tumor.currentLine.locus, tumor.contigTuple)
      for variant in list(tumorLocusVariantTable.keys()):
         if not matchedNormalVariantTable:  #we get false if the line was uncovered
            uncoveredNormalLine = UncoveredDataLine(tumor.currentLine.locus, variant, tumorLocusVariantTable[variant].variantType)
            outputLine = OutputLine(tumorLocusVariantTable[variant], uncoveredNormalLine)
            print(outputLine, file = outputFile)
         elif not variant in matchedNormalVariantTable:  #the tumor variant was not in the normal, but there was locus data returned (either a different variant or "None" if no variant was found in normal)
            normalDataGetter = list(matchedNormalVariantTable.keys())[0]
            unsupportedNormalLine = UnsupportedDataLine(matchedNormalVariantTable[normalDataGetter], variant, tumorLocusVariantTable[variant].variantType)
            outputLine = OutputLine(tumorLocusVariantTable[variant], unsupportedNormalLine)
            print(outputLine, file = outputFile)
         else:
            outputLine = OutputLine(tumorLocusVariantTable[variant], matchedNormalVariantTable[variant])
            print(outputLine, file = outputFile)
   tumorFile.close()
   normalFile.close()
   outputFile.close()
   if args.verbose:
      print("Compared %s loci in %s" %(locusCount, datetime.datetime.now() - startTime))