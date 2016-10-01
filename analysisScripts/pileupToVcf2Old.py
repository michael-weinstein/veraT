class CheckArgs(object):
   
   def __init__(self):
         import argparse
         import os
         import sys
         parser = argparse.ArgumentParser()
         parser.add_argument("-f", "--inputFile", help = "Input file")
         parser.add_argument("-o", "--outputFile", help = "Output file.")
         parser.add_argument("-n", "--minReadRequirement", help = "Require at least n reads supporting a variant for it to be emitted", type = int, default = 0)
         parser.add_argument("-p", "--minPercentageRequirement", help = "Require a variant to be supported by at least this percent of the reads", type = float, default = 0)
         parser.add_argument("-d", "--requireDoubleStranded", help = "Require a variant to be read on both strands to be emitted", action = 'store_true')
         parser.add_argument("-v", "--verbose", help = "Verbose mode", action = 'store_true')
         rawArgs = parser.parse_args()
         if rawArgs.inputFile:
            if os.path.isfile(rawArgs.inputFile):
                self.inputFile = rawArgs.inputFile
            else:
                raise RuntimeError("Input file not found: " + rawArgs.inputFile)
         else:
            raise RuntimeError("No input file specified.")
         self.outputFile = rawArgs.outputFile
         if self.outputFile == self.inputFile:
            raise RuntimeError("Error: Input and output files cannot be the same.")
         self.minReadRequirement = rawArgs.minReadRequirement
         self.requireDoubleStranded = rawArgs.requireDoubleStranded
         self.minPercentageRequirement = rawArgs.minPercentageRequirement
         if self.minPercentageRequirement >= 1 or self.minPercentageRequirement < 0:
            raise RuntimeError("Percentage requirement must less than 1 and greater than or equal to 0")
         self.verbose = rawArgs.verbose
   
class MPileupLine(object):
   
   def __init__(self, line, minReadRequirement = 0, minPercentageRequirement = 0, requireDoubleStranded = False, delimiter = "\t", autoparse = True):
      import re
      self.delimiter = delimiter
      self.minReadRequirement = minReadRequirement
      self.requireDoubleStranded = requireDoubleStranded
      self.minPercentageRequirement = minPercentageRequirement
      if self.minPercentageRequirement >= 1 or self.minPercentageRequirement < 0:
         raise RuntimeError("Percentage requirement must less than 1 and greater than or equal to 0")
      self.rawLine = line.strip()
      self.lineList = self.rawLine.split("\t")
      if len(self.lineList) < 4:
         self.isValidLine = False
      else:
         self.isValidLine = True
         self.contig = self.lineList[0]
         self.position = int(self.lineList[1])
         self.referenceBase = self.lineList[2].upper()
         self.readsCovering = int(self.lineList[3])
         self.isCovered = self.readsCovering > 0
         if self.isCovered:
            self.readBases = self.lineList[4]
            readMod = re.sub('\^.', "", self.readBases)
            self.readsCovering -= readMod.count("*")  #use the modified line so I can't have asterisks as map qualities in there, then remove remaining asterisks as they indicate previously counted gaps
            self.readQuality = self.lineList[5]
            self.isCovered = self.readsCovering > 0
            self.hasVariant = not set(readMod).issubset( {'.', ',', "$", "*"})
      if autoparse:
         if self.isValidLine and self.isCovered and self.hasVariant:
            self.parsePileupString()
         else:
            self.createNullOutputs()
      
   def parsePileupString(self):
      import sys
      counts = {}
      strands = {}
      position = 0
      varType = None
      self.variantCount = 0
      while position < len(self.readBases):
         character = self.readBases[position]
         if character == "*":
            position += 1
            continue
         if character in "$^":
            position += 1
            if character == "^":
               position += 1
            continue
         elif character in "atgcnATGCN.,":
            if self.referenceBase and not self.referenceBase in "Nn":
               if character in ".,":
                  if character == ".":
                     character = self.referenceBase.upper()
                  elif character == ",":
                     character = self.referenceBase.lower()
               else:
                  self.variantCount += 1
            variant = character.upper()
            position += 1
            if position < len(self.readBases):
               checkNext = self.readBases[position]
               if checkNext in "+-":
                  varType = checkNext
                  position += 1
                  numberEnd = position + 1
                  while numberEnd < len(self.readBases) and self.readBases[numberEnd].isdigit():
                     numberEnd += 1
                  indelSize = int(self.readBases[position : numberEnd])
                  position = numberEnd
                  indelSequence = self.readBases[position : position + indelSize]
                  position = position + indelSize
                  variant = (character.upper(), varType, indelSize, indelSequence.upper())
            if variant in counts:
               counts[variant] += 1
            else:
               counts[variant] = 1
               strands[variant] = ""  #initialize this now so we don't have to do it later.  Assumption is that any variant in counts should also be going into strands at the same time
            if character.isupper():
               strand = "+"
            elif character.islower():
               strand = "-"
            else:
               raise RuntimeError("Picked up an unhandled valid character %s analyzing line %s" %(character, self.line))
            strands[variant] += strand
         else:
            print("Skipping invalid character %s on line %s" %(character, self.line), file = sys.stdout)
      self.counts = counts
      self.strands = strands
      self.makeQualifiedVariantLines()
      
   def varToString(self, variant):
      if type(variant) == tuple:
         return "%s%d%s" %(variant[1], variant[2], variant[3])
      else:
         return variant
      
   def isIndel(self, variant, trueReturn = True, falseReturn = False):
      if len(variant) > 1:
         return trueReturn
      else:
         return falseReturn
   
   def makeQualifiedVariantLines(self):
      self.variantLineLists = []
      for key in list(self.counts.keys()):
         if type(key) == str and key.upper() == self.referenceBase.upper():
            continue
         readPercentage = self.counts[key]/self.readsCovering
         if not readPercentage >= self.minPercentageRequirement:
            continue
         if not self.counts[key] >= self.minReadRequirement:
            continue
         if self.requireDoubleStranded and not set(self.strands[key]) == set('+', '-'):
            continue
         indelFlag = self.isIndel(key, "INDEL;", "") + "NA=NA"
         thisVariantLineList = [self.counts[key], self.readsCovering, round(readPercentage, 4), self.contig, self.position, ".", self.referenceBase, self.varToString(key), "0", ".", indelFlag, "NA", "NA"]
         self.variantLineLists.append(thisVariantLineList)
      if not self.variantLineLists:
         self.variantLineStrings = []
         self.variantString = ""
         self.hasQualifiedVariants = False
      else:
         self.variantLineStrings = []
         for lineList in self.variantLineLists:
            stringed = [str(item) for item in lineList]
            self.variantLineStrings.append(self.delimiter.join(stringed))
         self.variantString = "\n".join(self.variantLineStrings)
         self.hasQualifiedVariants = True
         
   def createNullOutputs(self):
      self.variantLineLists = []
      self.variantLineStrings = []
      self.variantString = ""
      self.hasQualifiedVariants = False
      
   def __str__(self):
      return self.variantString
      


if __name__=="__main__":
   import datetime
   startTime = datetime.datetime.now()
   import sys
   args = CheckArgs()
   inputFileName = args.inputFile
   outputFileName = args.outputFile
   inputFile = open(inputFileName, 'r')
   if args.outputFile.upper == "STDOUT":
      usingStdout = True
      outputFile = sys.stdout
   else:
      usingStdout = False
      outputFile = open(outputFileName, 'w')
   inputLine = inputFile.readline().strip()
   progress = 0
   while inputLine:
      if args.verbose:
         if progress % 10000 == 0:
            print("Processed %s lines." %(progress), file = sys.stderr, end = "\r")
      progress += 1
      if inputLine.startswith("#"):
         inputLine = inputFile.readline().strip()
         continue
      line = MPileupLine(inputLine, args.minReadRequirement, args.minPercentageRequirement, args.requireDoubleStranded)
      if line.hasQualifiedVariants:
         print(line, file = outputFile)
      inputLine = inputFile.readline().strip()
   if not usingStdout:
      outputFile.close()
   inputFile.close()
   if args.verbose:
      print("Processed %s lines." %(progress), file = sys.stderr)
      print("Completed in %s" %(datetime.datetime.now() - startTime), file = sys.stderr)