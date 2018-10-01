#!/usr/bin/env python3

defaultLinesPerArrayJob = 2000000
killswitchFileName = "/u/home/a/alexhaol/killswitch"


class CheckArgs():  # class that checks arguments and ultimately returns a validated set of arguments to the main program

    def __init__(self):
        import argparse
        import os
        if killswitchFileName:
            if os.path.isfile(killswitchFileName):
                quit("Kill switch file found at %s" %killswitchFileName)
        parser = argparse.ArgumentParser()
        parser.add_argument("-f", "--mpileupFile", help="RNA mPileup", required=True)
        parser.add_argument("-s", "--somaticVariants", help="Pickle containing somatic variant analysis from DNA", required=True)
        parser.add_argument("-o", "--output", help="Output pickle file name")
        parser.add_argument("-m", "--minDiff", help="Minimum percent difference in expression vs. DNA mutant/wild-type ratios to consider worth scoring", default=10, type=int)
        parser.add_argument("-v", "--verbose", help="Verbose output mode", action='store_true')
        parser.add_argument("-p", "--noParallelChromosomes", help="Do not run chromosomes in parallel", action='store_true')
        parser.add_argument("--arrayJob", help="Indicates job is part of a job array.  Should only be entered by user for debugging purposes.", action='store_true')
        parser.add_argument("--mock", help="Do not actually submit jobs to queue", action='store_true')
        parser.add_argument("--noCleanup", help="Do not cleanup temporary directory when completed", action = 'store_true')
        parser.add_argument("--linesPerJob", help="Number of lines per array job", type = int, default = defaultLinesPerArrayJob)
        parser.add_argument("--workingDirectory", help="Array job working directory.")
        rawArgs = parser.parse_args()
        mpileupFile = rawArgs.mpileupFile
        if not os.path.isfile(mpileupFile):
            raise FileNotFoundError("Unable to find file %s" % mpileupFile)
        self.mpileupFile = os.path.abspath(mpileupFile)
        self.resumePickle = os.path.split(mpileupFile)[0] + os.sep + "mpileupRNAJobArray.pkl"
        self.resumeWaiting = os.path.isfile(self.resumePickle)
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
        self.parallelChromosomes = not rawArgs.noParallelChromosomes
        if rawArgs.arrayJob:
            workingDirectory = rawArgs.workingDirectory
            if not workingDirectory:
                raise ValueError("No working directory specified, but one is needed for array jobs to output properly.")
            self.arrayJob = int(os.environ["SGE_TASK_ID"])
            self.clockoutFile = workingDirectory + os.sep + "%s.done" %self.arrayJob
            self.output = workingDirectory + os.sep + "%s.pkl" %self.arrayJob
            self.parallelChromosomes = False
        else:
            self.clockoutFile = None
            self.arrayJob = False
        self.linesPerJob = rawArgs.linesPerJob
        self.mock = rawArgs.mock
        self.noCleanup = rawArgs.noCleanup
        if not self.arrayJob:
            self.skipLines = None
        else:
            self.skipLines = (self.arrayJob - 1) * self.linesPerJob


class ArrayJob(object):

    def __init__(self, jobNumber, workingDirectory):
        import os
        self.jobNumber = jobNumber
        self.outputFile = workingDirectory + os.sep + str(jobNumber) + ".pkl"
        self.clockoutFile = workingDirectory + os.sep + str(jobNumber) + ".done"
        self.completed = False

    def isDone(self):
        if self.completed:
            return True
        import os
        self.completed = os.path.isfile(self.clockoutFile)
        return self.completed


class ArrayJobs(object):

    def __init__(self, workingDirectory, args:CheckArgs):
        self.workingDirectory = workingDirectory
        self.args = args
        self.completed = False
        self.neededJobs = calculateNeededJobs(args.mpileupFile, args.linesPerJob)
        self.completedJobSet = set()
        self.remainingJobCount = self.neededJobs
        self.jobList = self.createJobList()

    def createJobList(self):
        jobList = []
        for jobNumber in range(1, self.neededJobs + 1):
            jobList.append(ArrayJob(jobNumber, self.workingDirectory))
        return jobList

    def isDone(self):
        allJobsDone = True
        for job in self.jobList:
            if job.isDone():
                self.completedJobSet.add(job.jobNumber)
                continue
            else:
                allJobsDone = False
        self.remainingJobCount = self.neededJobs - len(self.completedJobSet)
        self.completed = allJobsDone
        return self.completed

    def submit(self):
        import os
        import sys
        self.workingDirectory = os.path.abspath(self.workingDirectory)
        pythonInterpreter = sys.executable
        thisScript = os.path.abspath(__file__)
        arrayArgument = "-t 1-%s" %self.neededJobs
        mPileupCommand = [pythonInterpreter, thisScript, "--mpileupFile", self.args.mpileupFile, "--somaticVariants", self.args.somaticVariants, "--minDiff", self.args.minDiff, "--arrayJob", "--workingDirectory", self.workingDirectory]
        mPileupCommand = [str(item) for item in mPileupCommand]
        mPileupCommand = " ".join(mPileupCommand)
        qsubCommand = "qsub -cwd -V -N mpsubArray %s -l h_data=8G,time=2:00:00 -m a -o %s -e %s" % (arrayArgument, self.workingDirectory, self.workingDirectory)
        fullCommand = 'echo "%s" | %s' % (mPileupCommand, qsubCommand)
        submitted = False
        attempts = 0
        if self.args.mock:
            print("Mock submit:")
            print(fullCommand)
        else:
            while not submitted and attempts < 11:
                exitCode = os.system(fullCommand)
                submitted = exitCode == 0
                attempts += 1
            if not submitted:
                raise RuntimeError("Unable to submit job successfully")


def createTempDir(workingFolder, args:CheckArgs=False):  # makes a temporary directory for this run.  Completions will clock out here and results will be reported back to it.
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


def runScatterJobs(args:CheckArgs):
    import os
    import time
    import datetime
    outputDirectory = os.path.split(os.path.abspath(args.output))[0] + os.sep
    workingDirectory = createTempDir(outputDirectory, args)
    jobArray = ArrayJobs(workingDirectory, args)
    jobArray.submit()
    return waitOnJobCompletion(jobArray, args)


def waitOnJobCompletion(jobArray:ArrayJobs, args:CheckArgs):
    import time
    import datetime
    while not jobArray.isDone():
        startTime = datetime.datetime.now()
        if args.verbose:
            print("Awaiting %s of %s jobs       " % (jobArray.remainingJobCount, jobArray.neededJobs), end="\r")
        if (datetime.datetime.now() - startTime).seconds / 3600 >= 20:
            import pickle
            import sys
            resumePickleFileHandle = open(args.resumePickle, 'wb')
            pickle.dump(jobArray, resumePickleFileHandle)
            resumePickleFileHandle.close()
            print("Dumped job array pickle to %s. Exiting status 99" %args.resumePickle)
            sys.exit(99)
        time.sleep(300)
    if args.verbose:
        print("\nDONE!")
    return (jobArray, jobArray.workingDirectory)


def gatherResults(arrayJobs:ArrayJobs):
    import pickle
    rnaSupportTable = {}
    for job in arrayJobs.jobList:
        inputFile = open(job.outputFile, 'rb')
        partialData = pickle.load(inputFile)
        inputFile.close()
        for key in partialData:
            if partialData[key]["RNASupport"]:
                rnaSupportTable[key] = partialData[key]["RNASupport"]
    return rnaSupportTable


def openMPileupFile(filename):
    if filename.endswith(".gz"):
        import gzip
        return gzip.open(filename, 'rt')
    else:
        return open(filename, 'r')


def calculateNeededJobs(mpileupFile, linesPerJob):
    def countLinesInFile(fileHandle):
        lineNumber = 0
        for lineNumber, line in enumerate(fileHandle):
            pass
        return lineNumber + 1
    mpileup = openMPileupFile(mpileupFile)
    print("Counting lines in file")
    totalLines = countLinesInFile(mpileup)
    mpileup.close()
    neededJobs = -(-totalLines // linesPerJob)
    return neededJobs


def checkMPileupForRNASupport(mpileupFile, somaticVariants, somaticVariantTable, minDiff, verbose=False, arrayJob=False, skipLines=False, linesPerJob=None, directOutput = None):
    import variantDataHandler
    import scipy.stats
    import re
    contig = None #initializing this for display in progress reporter
    nonBaseRegex = re.compile("[^ATGC]", re.IGNORECASE)  # compiling this regex now so it only needs to be compiled once
    plusBaseRegex = re.compile("\+\d+[ATGC]+", re.IGNORECASE)
    mpileup = openMPileupFile(mpileupFile)
    supportData = {}
    lociOfInterest = [(item[0], str(item[1])) for item in somaticVariants]  # using a string of the position to improve performance, will convert hundreds of times now instead of converting millions of times during mpileup reading
    skipped = 0
    progress = 0
    startedRegionOfInterest = False
    #print("ChrRest: %s" % chromosomeRestriction)
    if skipLines:
        for i in range(skipLines):
            skipped += 1
            throwaway = mpileup.readline()
    for line in mpileup:
        if verbose:
            if progress % 10000 == 0:
                print("Processed %s lines. Currently on %s" %(progress, contig), end="\r")
            progress += 1
        if arrayJob:
            if progress > linesPerJob:
                break
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
            reads = re.sub(plusBaseRegex, "", reads)  # get rid of any plus (digit)(base) reads here, otherwise the base will remain as an alt
            reads = re.sub(nonBaseRegex, "", reads)  # get rid of anything not a base (we have already replaced the periods and commas)
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
                    supportData[foundHash] = variantDataHandler.RNASupportData(1, supportingDepthRNA, totalDepthRNA, None)
                    continue
                if totalDepthRNA < 10 or supportingDepthRNA < 3:
                    supportData[foundHash] = variantDataHandler.RNASupportData(2, supportingDepthRNA, totalDepthRNA, None)
                    continue
                if expressionDNARatio > (1 - minDiff / 100) and expressionDNARatio < (1 + minDiff / 100):  # anything not greater or less than minDiff percent off expected will be called as not significantly different
                    supportData[foundHash] = variantDataHandler.RNASupportData(4, supportingDepthRNA, totalDepthRNA, None)
                    continue
                oddsRatio, pvalue = scipy.stats.fisher_exact(
                    [[supportingDepthRNA, totalDepthRNA], [supportingDepthDNA, totalDepthDNA]])
                if pvalue > 0.05:
                    supportData[foundHash] = variantDataHandler.RNASupportData(4, supportingDepthRNA, totalDepthRNA, pvalue, oddsRatio)
                elif expressionDNARatio > 1:
                    supportData[foundHash] = variantDataHandler.RNASupportData(5, supportingDepthRNA, totalDepthRNA, pvalue, oddsRatio)
                else:
                    supportData[foundHash] = variantDataHandler.RNASupportData(3, supportingDepthRNA, totalDepthRNA, pvalue, oddsRatio)
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
    somaticsFile = open(args.somaticVariants, 'rb')
    somaticVariantTable = pickle.load(somaticsFile)
    somaticsFile.close()
    somaticVariants = list(somaticVariantTable.keys())
    for key in somaticVariants:
        somaticVariantTable[key]["RNASupport"] = variantDataHandler.RNASupportData(0, 0, 0)
    if args.parallelChromosomes:
        if not args.resumeWaiting:
            arrayJobs, tempdir = runScatterJobs(args)
        else:
            import pickle
            pickleFileHandle = open(args.resumePickle, 'rb')
            jobArray = pickle.load(pickleFileHandle)
            pickleFileHandle.close()
            arrayJobs, tempdir = waitOnJobCompletion(jobArray, args)
        rnaSupportTable = gatherResults(arrayJobs)
        if not args.noCleanup:
            import shutil
            print("All jobs completed, removing working directory at %s" %tempdir)
            shutil.rmtree(tempdir)
    else:
        rnaSupportTable = checkMPileupForRNASupport(args.mpileupFile, somaticVariants, somaticVariantTable, args.minDiff, args.verbose, args.arrayJob, args.skipLines, args.linesPerJob)
    for key in list(rnaSupportTable.keys()):
        somaticVariantTable[key]["RNASupport"] = rnaSupportTable[key]
    if args.output.upper().endswith(".PKL"):
        outputFile = open(args.output, 'wb')
        pickle.dump(somaticVariantTable, outputFile)
        outputFile.close()
    else:
        sortVariantDataTuples(somaticVariants)
        outputTable = createOutputTextTable(somaticVariants, somaticVariantTable)
        outputFile = open(args.output, 'w')
        for line in outputTable:
            print("\t".join(line), file=outputFile)
        outputFile.close()
    if args.clockoutFile:
        emptyFile = open(args.clockoutFile, 'w')
        emptyFile.close()
    quit()


if __name__ == '__main__':
    main()
