
import sys
import os
if not os.getcwd() in sys.path:
   sys.path.append(os.getcwd())

def readVariants(fname):
   variants = {}
   f = open(fname)
   for line in f:
      if line.startswith("#"): continue
      cols = line.strip().split("\t")
      if len(cols) < 4: continue
      chromNum = cols[3]
      loc = long(cols[4])
      variants[(chromNum, loc)] = line.strip()
   f.close()
   return variants

if __name__=="__main__":
   tumorVars = readVariants(sys.argv[1])
   normalVars = readVariants(sys.argv[2])
   missingVarString = "NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA"
   chromNums = set([chromNum for (chromNum, loc) in tumorVars])
   chromNums = list(chromNums)
   chromNums.sort()
   for chromNum in chromNums:
      locs = []
      for (c, loc) in tumorVars:
         if c == chromNum:
            locs.append(loc)
      locs.sort()
      for loc in locs:
         if normalVars.has_key((chromNum, loc)):
            print("%s\t%s" % (tumorVars[(chromNum, loc)], normalVars[(chromNum, loc)]))
         else:
            print("%s\t%s" % (tumorVars[(chromNum, loc)], missingVarString))
