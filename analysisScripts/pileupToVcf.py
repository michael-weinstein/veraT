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
            raise RunTimeError("No input file specified.")
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
   
   def __init__(self, line, minReadRequirement = 0, minPercentageRequirement = 0, requireDoubleStranded = False, delimiter = "\t"):
      import re
      self.delimiter = delimiter
      self.minReadRequirement = minReadRequirement
      self.requireDoubleStranded = requireDoubleStranded
      self.minPercentageRequirement = minPercentageRequirement
      if self.minPercentageRequirement >= 1 or self.minPercentageRequirement < 0:
         raise RuntimeError("Percentage requirement must less than 1 and greater than or equal to 0")
      self.rawLine = line.strip()
      self.lineList = rawLine.split("\t")
      if len(self.lineList) < 4:
         self.isValidLine = False
      else:
         self.isValidLine = True
         self.contig = lineList[0]
         self.position = int(lineList[1])
         self.referenceBase = lineList[2].upper()
         self.readsCovering = int(lineList[3])
         self.readBases = lineList[4]
         readMod = re.sub('\^.', "", self.readBases)
         self.readsCovering -= self.readMod.count("*")  #use the modified line so I can't have asterisks as map qualities in there, then remove remaining asterisks as they indicate previously counted gaps
         self.readQuality = lineList[5]
         self.isCovered = readsCovering > 0
         self.hasVariant = not set(readMod).issubset( {'.', ',', "$", "*"})
         self.variantCount = 0
      
   def parsePileupString(self):
      import sys
      counts = {}
      strands = {}
      position = 0
      varType = None
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
            checkNext = self.readBases[position]
            if checkNext in "+-":
               varType = checkNext
               position += 1
               numberEnd = position + 1
               while self.readBases[numberEnd].isdigit():
                  numberEnd += 1
               indelSize = int(self.readBases[position : numberEnd])
               position = numberEnd
               indelSequence = self.readBases[position : position + indelSize]
               position = position + indelSize
               variant = (character.upper(), varType, indelSize, indelSequence)
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
      
   def __str__(self):
      return self.variantString
      


if __name__=="__main__":
   import sys
   args = CheckArgs()
   inputFileName = args.inputFile
   outputFileName = args.outputFile
   inputFile = open(inputFileName, 'r')
   if args.outputFileName.upper == "STDOUT":
      usingStdout = True
      outputFile = sys.stdout
   else:
      usingStdout = False
      outputFile = open(outputFileName, 'w')
   inputLine = inputFile.readline()
   progress = 0
   while inputLine:
      if args.verbose:
         if progress % 10000 == 0:
            print("Processed %s lines.", file = sys.stderr, end = "\r")
      if inputLine.startswith("#"):
         inputLine = inputFile.readline()
         continue
      line = MPileupLine(inputLine, args.minReadRequirement, args.minPercentageRequirement, args.requireDoubleStranded)
      if line.isValidLine and line.isCovered and line.hasVariant:
         line.parsePileupString()
         if line.hasQualifiedVariants:
            print(line, file = outputFile)
      inputLine = inputFile.readLine()
   if not usingStdout:
      outputFile.close()
   inputFile.close()
         

      
'''
def updateStrands(curr, ch):
   forward = (ch == ch.upper())
   if curr == None:
      if forward: return "+"
      return "-"
   if curr == "-" and not forward:
      return curr
   if curr == "+" and forward:
      return curr
   return "+-"

def parsePileupString(st, refBase=None):
   # parse pileup string into counts of bases at this position
   pos = 0
   cts = {}
   strands = {}
   while pos < len(st):
      ch = st[pos]
      if refBase and refBase not in "nN" and ch in ".,":
         if ch == ".": ch = refBase.upper()
         else: ch = refBase.lower()
      if ch in "$^":
         pos += 1
         if ch == "^":
            pos += 1
      elif ch in "acgtnACGTN*":
         var = ch.upper()
         pos += 1
         if pos < len(st) and st[pos] in "-+":
            # insertion or deletion
            indelType = st[pos]
            pos += 1
            ll = 0
            while pos < len(st) and st[pos] in "0123456789":
               ll = ll*10 + int(st[pos])
               pos += 1
            indelChars = st[pos:(pos+ll)]
            pos += ll
            var = (ch.upper(), indelType, len(indelChars), indelChars.upper())
         cts[var] = cts.get(var, 0) + 1
         strands[var] = updateStrands(strands.get(var, None), ch)
      else:
         print("## Skipping character: " + ch + " at pos " + pos + " in string " + st, file = sys.stderr)
         pos += 1
   return (cts, strands)

def variantDesc(var):
   # either the new base or a string like '+2AC' for an indel
   if isinstance(var, tuple):
      return "%s%d%s" % (var[1], var[2], var[3])
   return str(var)

if __name__=="__main__":
   locs = set([])
   variantList = ["A", "C", "G", "T"]
   allResults = {}
   refBases = {}
   progress = 0
   for line in sys.stdin:
       if progress % 10000 == 0:
         print("Processed %s lines." %progress, file = sys.stderr, end = "\r")
       progress += 1
       cols = line.strip().split("\t")
       if len(cols) < 4:
          print("# Skipping line: " + line, file = sys.stderr)
          continue
       loc = (cols[0], int(cols[1]))
       locs.add(loc)
       refBase = cols[2].upper()
       if int(cols[3]) > 0:
          (cts, strands) = parsePileupString(cols[-2], refBase)
       else:
          (cts, strands) = ({}, {})
       for var in cts.keys():
          vd = variantDesc(var)
          if not vd in variantList: variantList.append(vd)
          allResults[(loc, vd)] = allResults.get((loc, vd), 0) + cts[var]
       if refBase != "N":
          refBases[loc] = refBase
       #acgt = map(lambda(x):cts.get(x,0), "ACGT")
       #totalCts = sum(cts.values())
       #refCts = cts.get(refBase, 0)
       # if totalCts - refCts >= 2:
       #print "\t".join(map(str, cols[:-1] + acgt + [sum(cts.values()) - sum(acgt), cts]))
   locList = list(locs)
   locList.sort()
   # print "Chrom\tPosition\tRefBase\t" + "\t".join(variantList)
   for loc in locList:
      totalCrossingReads = 0
      totalBases = 0
      refBase = refBases.get(loc, "N")
      for var in variantList:
         ct = allResults.get((loc, var), 0)
         totalCrossingReads += ct
         if var != "*": totalBases += ct
      for var in variantList:
         ct = allResults.get((loc, var), 0)
         if var != refBase and var != "*" and ct >= 2:
            indel=""
            if len(var)>1: indel="INDEL;"
            print("%d\t%d\t%.4f\t%s\t%s\t.\t%s\t%s\t0\t.\t%sNA=NA\tNA\tNA" % (ct, totalBases, 1.0*ct/totalBases, loc[0], loc[1], refBase, var, indel))
'''