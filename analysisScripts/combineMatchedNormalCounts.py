import sys
from pileupToVcf import updateStrands, parsePileupString

SUBSTITUTION = "Substitution"
INDEL = "Indel"

def parseInfoString(st):
   info = {}
   info["INDEL"] = SUBSTITUTION
   if st == "NA": return info
   fields = st.split(";")
   for field in fields:
      if field == "INDEL":
         info["INDEL"] = INDEL
         continue
      (name, val) = tuple(field.split("="))
      info[name] = val
   return info

def readVariantCounts(fname):
   chroms = {}
   f = open(fname)
   for line in f:
      if line.startswith("#"): continue
      cols = line.strip().split("\t")
      chromNum = cols[3]
      loc = long(cols[4])
      chroms[chromNum] = chroms.get(chromNum, 0) + 1
   f.close()
   return chroms

def readVariants(fname, theChromNum=None):
   variants = {}
   f = open(fname)
   for line in f:
      if line.startswith("#"): continue
      cols = line.strip().split("\t")
      chromNum = cols[3]
      if theChromNum and (chromNum != theChromNum): continue
      loc = long(cols[4])
      tumorVarInfo = parseInfoString(cols[10])
      normalVarInfo = parseInfoString(cols[23])
      tumorVarCts = cols[0]
      normalVarCts = cols[13]
      if tumorVarCts == "NA": tumorVarInfo["INDEL"] = "NA"
      if normalVarCts == "NA": normalVarInfo["INDEL"] = "NA"
      variants[(chromNum, loc)] = (tumorVarInfo.get("INDEL", False), normalVarInfo.get("INDEL", False), line.strip())
   f.close()
   return variants

def readPositionalCounts(fname, theChromNum=None):
   cts = {}
   f = open(fname)
   for line in f:
      cols = line.strip().split("\t")
      chromNum = cols[0]
      if theChromNum and (chromNum != theChromNum): continue
      loc = long(cols[1])
      if len(cols) == 4:
         line = line.strip() + "\t--\t--"
         (pileupDesc, pileupStrands) = ({}, {})
      else:
         (pileupDesc, pileupStrands) = parsePileupString(cols[4], cols[2].upper())
      cts[(chromNum, loc)] = (line.strip(), pileupDesc, pileupStrands)
   f.close()
   return cts

def processVariants(tumorAndNormalVars, tumorCts, matchedNormalCts, chromNum):
      missingVarString = "NA\tNA\tNA\t0\t--\t--"
      locs = []
      for (c, loc) in tumorAndNormalVars:
         if c == chromNum:
            locs.append(loc)
      locs.sort()
      for loc in locs:
         (tumorIsIndel, normalIsIndel, inf) = tumorAndNormalVars[(chromNum, loc)]
         (tumorInf, tumorPileupDesc, tumorPileupStrands) = tumorCts.get((chromNum, loc), (missingVarString, {}, {}))
         (normalInf, normalPileupDesc, normalPileupStrands) = matchedNormalCts.get((chromNum, loc), (missingVarString, {}, {}))
         # check pileup description?  get the variant at least.
         if tumorIsIndel == SUBSTITUTION:
            infCols = inf.strip().split("\t")
            oldBase = infCols[6]
            # use the most likely secondary base
            newBaseCount = -1
            newBase = None
            for k in "ACGT":
               if k != oldBase:
                  ct = tumorPileupDesc.get(k, 0)
                  if ct > newBaseCount:
                     newBaseCount = ct
                     newBase = k
            newBaseCount = tumorPileupDesc.get(newBase, 0)
            tumorCov = tumorInf.split("\t")[3]
            newNormalBaseCount = normalPileupDesc.get(newBase, 0)
            normalCov = normalInf.split("\t")[3]
            if newBaseCount > 0:
               chromName = str(infCols[3])
               if not chromName.startswith("chr"):
                  chromName = "chr" + chromName
               print("%s:%ld\t%s\t%s\t%s\t%d\t%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (chromName, long(infCols[4]), oldBase, newBase, "Substitution", newBaseCount, tumorCov, newNormalBaseCount, normalCov, tumorPileupStrands.get(newBase, "None"), normalPileupStrands.get(newBase, "None"), newBase, inf, tumorInf, normalInf))
         else:
            infCols = inf.strip().split("\t")
            # use the most likely indel
            newBaseCount = -1
            newBase = None
            for k in tumorPileupDesc.keys():
               if isinstance(k, tuple):
                  ct = tumorPileupDesc.get(k, 0)
                  if ct > newBaseCount:
                     newBaseCount = ct
                     newBase = k
            newBaseCount = tumorPileupDesc.get(newBase, 0)
            tumorCov = tumorInf.split("\t")[3]
            newNormalBaseCount = normalPileupDesc.get(newBase, 0)
            normalCov = normalInf.split("\t")[3]
            if newBaseCount > 0:
               oldBase = newBase[0]
               if newBase[1] == "+":
                  newBaseStr = "I-%d" % newBase[2]
               else:
                  newBaseStr = "D-%d" % newBase[2]
               newBaseChars = newBase[3]
               chromName = str(infCols[3])
               if not chromName.startswith("chr"):
                  chromName = "chr" + chromName
               print("%s:%ld\t%s\t%s\t%s\t%d\t%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (chromName, long(infCols[4]), oldBase, newBaseStr, "Indel", newBaseCount, tumorCov, newNormalBaseCount, normalCov, tumorPileupStrands.get(newBase, "None"), normalPileupStrands.get(newBase, "None"), newBaseChars, inf, tumorInf, normalInf))
            # print "%s\t%s\t%s\t%s\t%s" % (inf, tumorIsIndel, tumorInf, normalIsIndel, normalInf)

if __name__=="__main__":
   variantCounts = readVariantCounts(sys.argv[1])
   chromNums = set(variantCounts.keys())
   chromNums = list(chromNums)
   chromNums.sort()
   for chromNum in chromNums:
      print >> sys.stderr, "# Read ", variantCounts[chromNum], " variants on chrom ", chromNum, " from ", sys.argv[1]
   for chromNum in chromNums:
      tumorAndNormalVars = readVariants(sys.argv[1], chromNum)
      tumorCts = readPositionalCounts(sys.argv[2], chromNum)
      matchedNormalCts = readPositionalCounts(sys.argv[3], chromNum)      
      processVariants(tumorAndNormalVars, tumorCts, matchedNormalCts, chromNum)
