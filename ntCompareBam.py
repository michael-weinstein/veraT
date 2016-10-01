#!/usr/bin/env python3

class CheckArgs(object):
    
    def __init__(self):
        import argparse
        import os
        inputFileArgs = {"normal" : "Normal paired-end 1",
                         "tumor" : "Tumor paired-end 1"}
        inputFileList = list(inputFileArgs.keys())
        parser = argparse.ArgumentParser()
        for flag in inputFileList:
            parser.add_argument("--" + flag, help = inputFileArgs[flag])
        parser.add_argument("-d", "--outputDirectory", help = "Specify an output directory", default = "")
        parser.add_argument("-m", "--mock", help = "Mock submission (don't actually submit to scheduler)", action = 'store_true')
        parser.add_argument("-9", "--clobber", help = "Overwrite existing files without asking.", action = 'store_true')
        parser.add_argument("-n", "--sampleName", help = "The unique sample identifier. Seriously, make it unique.  Or not... be lazy and reap the consequences.")
        parser.add_argument("-r", "--referenceGenome", help = "The reference genome fasta")
        parser.add_argument("-M", "--emailAddress", help = "Email address for notifications", default = "")
        parser.add_argument("--emailConditions", help = "Conditions to email for (b)egin, (e)nd, (a)bort.", default = "a")
        parser.add_argument("--coverageRequirement", help = "Required coverage to emit a variant from tumor or normal (matched sites from normal will bypass this requirement).", type = int, default = 2)
        rawArgs = parser.parse_args()
        self.mock = rawArgs.mock
        clobber = rawArgs.clobber
        if not clobber:
            self.clobber = None  #using none as a default type here.  We will force the user to set this if it becomes an issue.  Hopefully it won't because they are using unique names.
        else:
            self.clobber = True
        inputFilePaths = {}
        inputFilePathList = []   #using this to check for duplicates.  All files should be unique.
        for flag in inputFileList:
            currentFilePath = rawArgs.__dict__[flag]
            if not currentFilePath:
                if flag.endswith("2"):
                    self.__dict__[flag] = ""
                else:
                    raise RuntimeError("No %s argument passed, but a %s is required." %(flag, inputFileArgs[flag]))
            else:
                if currentFilePath in inputFilePathList:
                    raise RuntimeError("File path %s was given as two different files.  All files should be unique." %(currentFilePath))
                inputFilePathList.append(currentFilePath)
                if not os.path.isfile(currentFilePath):
                    raise RuntimeError("File path for %s was given as %s, but the file could not be found." %(inputFileArgs[flag].lower(), currentFilePath))
                self.__dict__[flag] = currentFilePath
        if not rawArgs.sampleName:
            raise RuntimeError("No sample name was specified under the -n or --sampleName arguments.")
        else:
            self.sampleName = rawArgs.sampleName
        referenceGenome = rawArgs.referenceGenome
        if not referenceGenome:
            raise RuntimeError("You must specify a reference genome")
        if not os.path.isfile(referenceGenome):
            raise FileNotFoundError("Reference genome fasta file was not found. Expected file: %s" %(referenceGenome))
        self.referenceGenome = referenceGenome
        self.emailAddress = rawArgs.emailAddress
        self.emailConditions = rawArgs.emailConditions
        self.outputDirectory = rawArgs.outputDirectory
        self.coverageRequirement = rawArgs.coverageRequirement
            
class JobList():
    def __init__(self):
        from runners import programRunners
        from runners import houseKeeping
        #deduplicate all the things
        self.deduplicateNormal = programRunners.Deduplicate(args.sampleName + ".normal", args.normal, clobber = args.clobber, outputDirectory = args.outputDirectory)
        self.deduplicateTumor = programRunners.Deduplicate(args.sampleName + ".tumor", args.tumor, clobber = self.deduplicateNormal.clobber, outputDirectory = args.outputDirectory)
        #mPileup all the things
        self.mPileupNormal = programRunners.MPileup(args.sampleName + ".normal", self.deduplicateNormal.bamOut, args.referenceGenome, clobber = self.deduplicateTumor.clobber, outputDirectory = args.outputDirectory)
        self.mPileupTumor = programRunners.MPileup(args.sampleName + ".tumor", self.deduplicateTumor.bamOut, args.referenceGenome, clobber = self.mPileupNormal.clobber, outputDirectory = args.outputDirectory)
        #extract variants from all the tumors
        self.extractTumor = programRunners.ExtractVariantsTumor(args.sampleName + ".tumor", self.mPileupTumor.mPileupOut, minSupport = args.coverageRequirement, requireDoubleStranded = False, outputDirectory = args.outputDirectory, clobber = self.mPileupTumor.clobber)
        #delete mPileups
        self.deleteNormalPileup = houseKeeping.Delete(self.mPileupNormal.mPileupOut)
        self.deleteTumorPileup = houseKeeping.Delete(self.mPileupTumor.mPileupOut)
        #extract matched variant lines from the normal for comparisons
        self.extractNormalVsTumor = programRunners.ExtractVariantsNormal(args.sampleName + ".norm2tumor", self.mPileupNormal.mPileupOut, self.extractTumor.targetList, "normal2tumor", minSupport = args.coverageRequirement, clobber = self.extractTumor.clobber, outputDirectory = args.outputDirectory)
        #generate combined outputs
        self.combineNormalVsTumor = programRunners.CombineExtractedVariants(args.sampleName + ".tumorVcf", self.extractTumor.variantsOut, self.extractNormalVsTumor.variantsOut , "tumorToNormal", self.extractNormalVsTumor.clobber, args.outputDirectory)
        self.deleteDedupBamNormal = houseKeeping.Delete(self.deduplicateNormal.bamOut)
        self.deleteDedupBamTumor = houseKeeping.Delete(self.deduplicateTumor.bamOut)
        self.commandList = [self.deduplicateNormal.deduplicateCommand,
                            self.deduplicateTumor.deduplicateCommand,
                            self.mPileupNormal.mPileupCommand,
                            self.mPileupTumor.mPileupCommand,
                            self.extractTumor.extractCommand,
                            self.extractNormalVsTumor.extractCommand,
                            self.combineNormalVsTumor.combineCommand,
                            self.deleteNormalPileup.deleteCommand,
                            self.deleteTumorPileup.deleteCommand,
                            self.deleteDedupBamNormal.deleteCommand,
                            self.deleteDedupBamTumor.deleteCommand]
        
def submitJobs(jobList):
    from runners import genericRunners
    tempDir = genericRunners.createTempDir(args.sampleName)
    writeCommandList(jobList.commandList, tempDir)
    deduplicateNormaldeduplicate = genericRunners.HoffmanJob(1, [], args.sampleName + "dn", tempDir, args.emailAddress, args.emailConditions, mock = args.mock)
    deduplicateTumordeduplicate = genericRunners.HoffmanJob(2, [], args.sampleName + "dt", tempDir, args.emailAddress, args.emailConditions, mock = args.mock)
    mPileupNormalmPileup = genericRunners.HoffmanJob(3, [deduplicateNormaldeduplicate.jobID], args.sampleName + "mn", tempDir, args.emailAddress, args.emailConditions, mock = args.mock)
    mPileupTumormPileup = genericRunners.HoffmanJob(4, [deduplicateTumordeduplicate.jobID], args.sampleName + "mt", tempDir, args.emailAddress, args.emailConditions, mock = args.mock)
    extractTumorextract = genericRunners.HoffmanJob(5, [mPileupTumormPileup.jobID], args.sampleName + "et", tempDir, args.emailAddress, args.emailConditions, mock = args.mock)
    extractNormalVsTumorextract = genericRunners.HoffmanJob(6, [mPileupNormalmPileup.jobID, extractTumorextract.jobID], args.sampleName + "nvt", tempDir, args.emailAddress, args.emailConditions, mock = args.mock)
    combineTumorVsNormal = genericRunners.HoffmanJob(7, [extractTumorextract.jobID, extractNormalVsTumorextract.jobID], args.sampleName + "tvcf", tempDir, args.emailAddress, "ea", mock = args.mock)
    deleteNormalPileup = genericRunners.HoffmanJob(8, [extractNormalVsTumorextract.jobID], args.sampleName + "delnmp", tempDir, args.emailAddress, args.emailConditions, mock = args.mock)
    deleteTumorPileup = genericRunners.HoffmanJob(9, [extractTumorextract.jobID], args.sampleName + "deltmp", tempDir, args.emailAddress, args.emailConditions, mock = args.mock)
    deleteDedupBamNormal = genericRunners.HoffmanJob(10, [mPileupNormalmPileup.jobID], args.sampleName + "delndb", tempDir, args.emailAddress, args.emailConditions, mock = args.mock)
    deleteDedupBamTumor = genericRunners.HoffmanJob(11, [mPileupTumormPileup.jobID], args.sampleName + "deltdb", tempDir, args.emailAddress, args.emailConditions, mock = args.mock)
    
def writeCommandList(commandList, tempDir):
    import os
    import pickle
    commandListFileName = tempDir + os.sep + "jobs.pkl"
    commandListFile = open(commandListFileName, 'wb')
    commandListForDump = ["Zeroth job place holder"] #hoffman indexes array jobs to 1, so there can be no zeroth job
    commandListForDump += commandList
    pickle.dump(commandListForDump, commandListFile)
    commandListFile.close()
    
def main():
    global args
    args = CheckArgs()
    if args.outputDirectory:
        import os
        try:
            os.mkdir(args.outputDirectory)
        except OSError:
            pass
    jobList = JobList()
    submitJobs(jobList)
    print("Completed")
    
main()