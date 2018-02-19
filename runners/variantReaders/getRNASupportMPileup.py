#!/usr/bin/env python3

class CheckArgs():  # class that checks arguments and ultimately returns a validated set of arguments to the main program

    def __init__(self):
        import argparse
        import os
        parser = argparse.ArgumentParser()
        parser.add_argument("-f", "--mpileupFile", help="RNA mPileup", required=True)
        parser.add_argument("-s", "--somaticVariants", help="Pickle containing somatic variant analysis from DNA",
                            required=True)
        parser.add_argument("-o", "--output", help="Output pickle file name")
        parser.add_argument("-m", "--minDiff",
                            help="Minimum percent difference in expression vs. DNA mutant/wild-type ratios to consider worth scoring",
                            default=10, type=int)
        parser.add_argument("-v", "--verbose", help="Verbose output mode", action='store_true')
        parser.add_argument("-p", "--parallelChromosomes", help="Run chromosomes in parallel", action='store_true')
        parser.add_argument("--chromosome",
                            help="Specifies the chromosome for analysis.  Should only be used manually for debug.")
        parser.add_argument("--clockoutFile", help="Clockout file for scatter/gather jobs")
        parser.add_argument("--mock", help="Do not actually submit jobs to queue", action='store_true')
        rawArgs = parser.parse_args()
        mpileupFile = rawArgs.mpileupFile
        if not os.path.isfile(mpileupFile):
            raise FileNotFoundError("Unable to find file %s" % mpileupFile)
        self.mpileupFile = mpileupFile
        somaticVariants = rawArgs.somaticVariants
        if not os.path.isfile(somaticVariants):
            raise FileNotFoundError("Unable to find file %s" % somaticVariants)
        self.somaticVariants = somaticVariants
        output = rawArgs.output
        if not output:
            output = self.mpileupFile + ".rnaSupport.pkl"
        self.output = output
        self.minDiff = rawArgs.minDiff
        self.verbose = rawArgs.verbose
        self.parallelChromosomes = rawArgs.parallelChromosomes
        self.chromosome = rawArgs.chromosome
        self.clockoutFile = rawArgs.clockoutFile
        self.mock = rawArgs.mock


class ChromosomeIndex(object):

    def __init__(self, chromosome, index):
        self.chromosome = chromosome
        self.index = index
        self.argument = "%s,%s" % (self.chromosome, self.index)

    def __str__(self):
        return self.argument


class ScatterJob(object):

    def __init__(self, chromosome, workingDirectory, args):
        import os
        self.chromosome = chromosome
        self.workingDirectory = workingDirectory
        self.outputFile = workingDirectory + os.sep + chromosome.chromosome + ".pkl"
        self.clockoutFile = workingDirectory + os.sep + chromosome.chromosome + ".done"
        self.args = args
        self.completed = False

    def submit(self, mock=False):
        import os
        import sys
        pythonInterpreter = sys.executable
        thisScript = os.path.abspath(__file__)
        mPileupCommand = [pythonInterpreter, thisScript, "--mpileupFile", args.mpileupFile, "--somaticVariants",
                          args.somaticVariants, "--output", self.outputFile, "--minDiff", args.minDiff, "--chromosome",
                          self.chromosome.argument, "--clockoutFile", self.clockoutFile]
        qsubCommand = "qsub -cwd -V -N mpsub%s -l h_data=8G,time=24:00:00 -m a" % self.chromosome.chromosome
        fullCommand = 'echo "%s" | %s' % (mPileupCommand, qsubCommand)
        submitted = False
        attempts = 0
        if mock:
            print("Mock submit:")
            print(fullCommand)
        else:
            while not submitted and attempts < 11:
                exitCode = os.system(fullCommand)
                submitted = exitCode == 0
                attempts += 1
            if not submitted:
                raise RuntimeError("Unable to submit job successfully")

    def checkCompletion(self):
        import os
        if os.path.isfile(self.clockoutFile):
            self.completed = True
            return True
        else:
            return False


def createTempDir(workingFolder,
                  args=False):  # makes a temporary directory for this run.  Completions will clock out here and results will be reported back to it.
    if args and args.verbose:  # If the user wants reports on what we're doing...
        print("Creating temporary directory")  # tell them
    import re  # import the library for regular expressions
    import os  # import the library for OS system calls
    import datetime  # import the library for date and time usage
    successful = False  # initialize a value
    while not successful:  # until the value of success is true
        currenttime = datetime.datetime.now()  # get the current date and time
        currenttime = str(currenttime)  # convert the datetime object to a string
        currenttime = re.sub(r'\W', '',
                             currenttime)  # get rid of any characters that are not numbers, letters, or underscores by replacing them with nothing
        tempdir = workingFolder + os.sep + '.RNAmpAnalysis' + currenttime  # create a string with the name of my temporary directory
        if os.path.isdir(tempdir):  # if there's already a directory with the name of our desired tempdir
            continue  # go back and try again
        try:
            os.mkdir(tempdir)  # try making the temporary directory
        except OSError:  # if it doesn't work
            continue  # go back and try again
        successful = True  # otherwise, set successful to true, return to the start of the loop, and exit
    return tempdir


def runScatterJobs(args):
    import os
    import time
    chromosomeIndex = createChromosomeIndex(args.mpileupFile)
    scatterJobs = []
    outputDirectory = os.path.split(os.path.abspath(args.output))[0] + os.sep
    workingDirectory = createTempDir(outputDirectory, args)
    for chromosome in chromosomeIndex:
        scatterJobs.append(ScatterJob(chromosome, workingDirectory, args))
    for job in scatterJobs:
        job.submit(args.mock)
    completed = False
    while not completed:
        completed = True
        for job in scatterJobs:
            if not job.completed:
                if not job.checkCompletion():
                    completed = False
        if args.verbose:
            completedJobs = 0
            for job in scatterJobs:
                if job.completed:
                    completedJobs += 1
            print("Awaiting %s of %s jobs       " % (len(scatterJobs) - completedJobs, len(scatterJobs)), end="\r")
        time.sleep(5)
    if args.verbose:
        print("\nDONE!")
    return scatterJobs


def gatherResults(scatterJobs):
    import pickle
    rnaSupportTable = {}
    for job in scatterJobs:
        inputFile = open(job.outputFile, 'rb')
        partialData = pickle.load(inputFile)
        inputFile.close()
        rnaSupportTable.update(partialData)
    return rnaSupportTable


def openMPileupFile(filename):
    if filename.endswith(".gz"):
        import gzip
        return gzip.open(filename, 'rt')
    else:
        return open(filename, 'r')


def createChromosomeIndex(mpileupFile):
    import re
    currentLine = 0
    chromosomeList = []
    mpileup = openMPileupFile(mpileupFile)
    lastContig = None
    contigRegex = re.compile("^(.+?)\t")
    for line in mpileup:
        chrGrab = re.match(contigRegex, line)
        chromosome = chrGrab.group(1)
        if not chromosome == lastContig:
            lastContig = chromosome
            chromosomeList.append(ChromosomeIndex(chromosome, currentLine))
        currentLine += 1
    mpileup.close()
    return chromosomeList


def checkMPileupForRNASupport(mpileupFile, somaticVariants, somaticVariantTable, minDiff, verbose=False,
                              chromosomeRestriction=False):
    import variantDataHandler
    import scipy.stats
    import re
    nonBaseRegex = re.compile("[^ATGC]", re.IGNORECASE)  # compiling this regex now so it only needs to be compiled once
    plusBaseRegex = re.compile("\+\d+[ATGC]+", re.IGNORECASE)
    contigRegex = re.compile("^(.+?)\t")
    mpileup = openMPileupFile(mpileupFile)
    supportData = {}
    lociOfInterest = [(item[0], str(item[1])) for item in
                      somaticVariants]  # using a string of the position to improve performance, will convert hundreds of times now instead of converting millions of times during mpileup reading
    progress = 0
    startedRegionOfInterest = False
    if chromosomeRestriction and "," in chromosomeRestriction:
        chromosomeRestriction, jumpLines = chromosomeRestriction.split(",")
        for i in range(jumpLines):
            throwaway = mpileup.readline()
    for line in mpileup:
        if chromosomeRestriction:
            chrGrab = re.match(contigRegex, line)
            chromosome = chrGrab.group(1)
            inRegionOfInterest = chromosome == chromosomeRestriction
            if not inRegionOfInterest:
                if not started:
                    continue
                else:
                    break
            else:
                started = True
        if verbose:
            if progress % 10000 == 0:
                print("Processed %s lines" % progress, end="\r")
            progress += 1
        line = line.strip()
        if not line:
            continue
        elif line.startswith("#"):
            continue
        else:
            try:
                contig, position, ref, depth, reads, quality = line.split()
            except ValueError:  # when lines with no reads are included, we get 4 columns (reads and qualities are missing)
                continue
            locus = (contig, position)
            if not locus in lociOfInterest:
                continue  # if we read past this, we know we are looking at a locus with a somatic change
            position = int(position)
            totalDepthRNA = int(depth)
            reads = reads.upper().replace(".", ref).replace(",", ref)
            reads = re.sub(plusBaseRegex, "",
                           reads)  # get rid of any plus (digit)(base) reads here, otherwise the base will remain as an alt
            reads = re.sub(nonBaseRegex, "",
                           reads)  # get rid of anything not a base (we have already replaced the periods and commas)
            foundHashesAtSite = []
            for variant in somaticVariants:
                if contig == variant[0] and position == variant[1]:
                    foundHashesAtSite.append(variant)
            for foundHash in foundHashesAtSite:
                if not (len(foundHash[2]) == 1 and len(foundHash[3]) == 1):  # indel catcher
                    continue
                altAllele = foundHash[3]
                supportingDepthRNA = reads.count(altAllele)
                totalDepthDNA = somaticVariantTable[foundHash]["combined"].tumorDepth
                supportingDepthDNA = somaticVariantTable[foundHash]["combined"].tumorSupporting
                expressionDNARatio = (supportingDepthRNA / totalDepthRNA) / (supportingDepthDNA / totalDepthDNA)
                if not supportingDepthRNA:
                    supportData[foundHash] = variantDataHandler.RNASupportData(1, supportingDepthRNA, totalDepthRNA,
                                                                               None)
                    continue
                if totalDepthRNA < 10 or supportingDepthRNA < 3:
                    supportData[foundHash] = variantDataHandler.RNASupportData(2, supportingDepthRNA, totalDepthRNA,
                                                                               None)
                    continue
                if expressionDNARatio > (1 - minDiff / 100) and expressionDNARatio < (
                        1 + minDiff / 100):  # anything not greater or less than minDiff percent off expected will be called as not significantly different
                    supportData[foundHash] = variantDataHandler.RNASupportData(4, supportingDepthRNA, totalDepthRNA,
                                                                               None)
                    continue
                oddsRatio, pvalue = scipy.stats.fisher_exact(
                    [[supportingDepthRNA, totalDepthRNA], [supportingDepthDNA, totalDepthDNA]])
                if pvalue > 0.05:
                    supportData[foundHash] = variantDataHandler.RNASupportData(4, supportingDepthRNA, totalDepthRNA,
                                                                               pvalue, oddsRatio)
                elif expressionDNARatio > 1:
                    supportData[foundHash] = variantDataHandler.RNASupportData(5, supportingDepthRNA, totalDepthRNA,
                                                                               pvalue, oddsRatio)
                else:
                    supportData[foundHash] = variantDataHandler.RNASupportData(3, supportingDepthRNA, totalDepthRNA,
                                                                               pvalue, oddsRatio)
    if verbose:
        print("Processed %s lines" % progress)
    mpileup.close()
    return supportData


def createOutputTextTable(sortedAcceptedVariantInfoTuples, variantDicts):
    sources = sorted(list(variantDicts[sortedAcceptedVariantInfoTuples[0]].keys()))
    sources.remove("hits")
    sources.remove("combined")
    sources.remove("RNASupport")
    headerLine = ["contig", "position", "ref_allele", "alt_allele", "CountMethod"] + sources + ["RNASupport"]
    outputTable = [headerLine]
    for variant in sortedAcceptedVariantInfoTuples:
        contig, position, ref, alt = variant
        outputLine = [contig, position, ref, alt]
        countData = []
        hits = variantDicts[variant]["hits"]
        for source in sources:
            if source in variantDicts[variant] and variantDicts[variant][source]:
                countData.append(variantDicts[variant][source].readCountString)
            else:
                countData.append("NA")
        countData.append(str(variantDicts[variant]["RNASupport"]))
        outputLine = outputLine + [hits] + countData
        outputLine = (str(field) for field in outputLine)
        outputTable.append(outputLine)
    return outputTable


def sortVariantDataTuples(variantDataTupleList):
    import operator
    import variantSupport
    # goal is to have variants sorted first by contig in the standard order, then by position within the contig.  This will be done by first sorting on position and then on contig using a function that specifies the order.
    variantDataTupleList.sort(key=operator.itemgetter(1))
    variantDataTupleList.sort(key=variantSupport.contigSortValueFromSomaticTuple)


def main():
    args = CheckArgs()
    import pickle
    import variantDataHandler
    if args.parallelChromosomes:
        scatterJobs = runScatterJobs(args)
        rnaSupportTable = gatherResults(scatterJobs)
    else:
        somaticsFile = open(args.somaticVariants, 'rb')
        somaticVariantTable = pickle.load(somaticsFile)
        somaticsFile.close()
        somaticVariants = list(somaticVariantTable.keys())
        for key in somaticVariants:
            somaticVariantTable[key]["RNASupport"] = variantDataHandler.RNASupportData(0, 0, 0)
        rnaSupportTable = checkMPileupForRNASupport(args.mpileupFile, somaticVariants, somaticVariantTable,
                                                    args.minDiff, args.verbose)
    for key in list(rnaSupportTable.keys()):
        somaticVariantTable[key]["RNASupport"] = rnaSupportTable[key]
    sortVariantDataTuples(somaticVariants)
    if args.output.upper().endswith(".PKL"):
        outputFile = open(args.output, 'wb')
        pickle.dump(somaticVariantTable, outputFile)
        outputFile.close()
    else:
        outputTable = createOutputTextTable(somaticVariants, somaticVariantTable)
        outputFile = open(args.output, 'w')
        for line in outputTable:
            print("\t".join(line), file=outputFile)
        outputFile.close()
    quit()


if __name__ == '__main__':
    main()
