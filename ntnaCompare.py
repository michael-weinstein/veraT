#!/usr/bin/env python3

class CheckArgs(object):
    
    def __init__(self):
        import argparse
        import os
        inputFileArgs = {"normal1" : "Normal paired-end 1",
                         "normal2" : "Normal paired-end 2",
                         "tumor1" : "Tumor paired-end 1",
                         "tumor2" : "Tumor paired-end 2",
                         "avatar1" : "Avatar paired-end 1",
                         "avatar2" : "Avatar paired-end 2",
                         "neurosphere1" : "Neurosphere paired-end 1",
                         "neurosphere2" : "Neurosphere paired-end 2"}
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
        parser.add_argument("--bwaCores", help = "Cores to give BWA", type = int, default = 2)
        parser.add_argument("--coverageRequirement", help = "Required coverage to emit a variant from tumor or normal (matched sites from normal will bypass this requirement).", type = int, default = 0)
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
        self.bwaCores = rawArgs.bwaCores
        self.coverageRequirement = rawArgs.coverageRequirement
            
class JobList():
    def __init__(self):
        from runners import programRunners
        from runners import houseKeeping
        #align all the things
        self.alignNormal = programRunners.BWAlign(args.sampleName + ".normal", args.referenceGenome, args.normal1, args.normal2, cores = args.bwaCores, clobber = args.clobber, outputDirectory = args.outputDirectory)
        self.alignTumor = programRunners.BWAlign(args.sampleName + ".tumor", args.referenceGenome, args.tumor1, args.tumor2, cores = args.bwaCores, clobber = self.alignNormal.clobber, outputDirectory = args.outputDirectory)
        self.alignNeurosphere = programRunners.BWAlign(args.sampleName + ".neurosphere", args.referenceGenome, args.neurosphere1, args.neurosphere2, cores = args.bwaCores, clobber = self.alignTumor.clobber, outputDirectory = args.outputDirectory)
        self.alignAvatar = programRunners.BWAlign(args.sampleName + ".avatar", args.referenceGenome, args.avatar1, args.avatar2, cores = args.bwaCores, clobber = self.alignNeurosphere.clobber, outputDirectory = args.outputDirectory)
        #delete sai files
        self.deleteNormalSai = houseKeeping.Delete([self.alignNormal.pe1Out, self.alignNormal.pe2Out])
        self.deleteTumorSai = houseKeeping.Delete([self.alignTumor.pe1Out, self.alignTumor.pe2Out])
        self.deleteNeurosphereSai = houseKeeping.Delete([self.alignNeurosphere.pe1Out, self.alignNeurosphere.pe2Out])
        self.deleteAvatarSai = houseKeeping.Delete([self.alignAvatar.pe1Out, self.alignAvatar.pe2Out])
        #compress all the things
        self.makeBAMNormal = programRunners.SAMtoBAM(args.sampleName + ".normal", self.alignNormal.samOut, clobber = self.alignAvatar.clobber, outputDirectory = args.outputDirectory)  #clobber can be changed during the previous object initialization.  Using that to keep current on changes.
        self.makeBAMTumor = programRunners.SAMtoBAM(args.sampleName + ".tumor", self.alignTumor.samOut, clobber = self.makeBAMNormal.clobber, outputDirectory = args.outputDirectory)  #clobber can be changed during the previous object initialization.  Using that to keep current on changes.
        self.makeBAMNeurosphere = programRunners.SAMtoBAM(args.sampleName + ".neurosphere", self.alignNeurosphere.samOut, clobber = self.makeBAMTumor.clobber, outputDirectory = args.outputDirectory)  #clobber can be changed during the previous object initialization.  Using that to keep current on changes.
        self.makeBAMAvatar = programRunners.SAMtoBAM(args.sampleName + ".avatar", self.alignAvatar.samOut, clobber = self.makeBAMNeurosphere.clobber, outputDirectory = args.outputDirectory)  #clobber can be changed during the previous object initialization.  Using that to keep current on changes.
        #delete sam files
        self.deleteNormalSam = houseKeeping.Delete(self.alignNormal.samOut)
        self.deleteTumorSam = houseKeeping.Delete(self.alignTumor.samOut)
        self.deleteNeurosphereSam = houseKeeping.Delete(self.alignNeurosphere.samOut)
        self.deleteAvatarSam = houseKeeping.Delete(self.alignAvatar.samOut)
        #deduplicate all the things
        self.deduplicateNormal = programRunners.Deduplicate(args.sampleName + ".normal", self.makeBAMNormal.bamOut, clobber = self.makeBAMAvatar.clobber, outputDirectory = args.outputDirectory)
        self.deduplicateTumor = programRunners.Deduplicate(args.sampleName + ".tumor", self.makeBAMTumor.bamOut, clobber = self.deduplicateNormal.clobber, outputDirectory = args.outputDirectory)
        self.deduplicateNeurosphere = programRunners.Deduplicate(args.sampleName + ".neurosphere", self.makeBAMNeurosphere.bamOut, clobber = self.deduplicateTumor.clobber, outputDirectory = args.outputDirectory)
        self.deduplicateAvatar = programRunners.Deduplicate(args.sampleName + ".avatar", self.makeBAMAvatar.bamOut, clobber = self.deduplicateNeurosphere.clobber, outputDirectory = args.outputDirectory)
        #delete non-deuplicated BAMs
        self.deleteFirstNormalBam = houseKeeping.Delete(self.makeBAMNormal.bamOut)
        self.deleteFirstTumorBam = houseKeeping.Delete(self.makeBAMTumor.bamOut)
        self.deleteFirstNeurosphereBam = houseKeeping.Delete(self.makeBAMNeurosphere.bamOut)
        self.deleteFirstAvatarBam = houseKeeping.Delete(self.makeBAMAvatar.bamOut)
        #mPileup all the things
        self.mPileupNormal = programRunners.MPileup(args.sampleName + ".normal", self.deduplicateNormal.bamOut, args.referenceGenome, clobber = self.deduplicateAvatar.clobber, outputDirectory = args.outputDirectory)
        self.mPileupTumor = programRunners.MPileup(args.sampleName + ".tumor", self.deduplicateTumor.bamOut, args.referenceGenome, clobber = self.mPileupNormal.clobber, outputDirectory = args.outputDirectory)
        self.mPileupNeurosphere = programRunners.MPileup(args.sampleName + ".neurosphere", self.deduplicateNeurosphere.bamOut, args.referenceGenome, clobber = self.mPileupTumor.clobber, outputDirectory = args.outputDirectory)
        self.mPileupAvatar = programRunners.MPileup(args.sampleName + ".avatar", self.deduplicateAvatar.bamOut, args.referenceGenome, clobber = self.mPileupNeurosphere.clobber, outputDirectory = args.outputDirectory)
        #extract variants from all the tumors
        self.extractTumor = programRunners.ExtractVariantsTumor(args.sampleName + ".tumor", self.mPileupTumor.mPileupOut, minSupport = args.coverageRequirement, requireDoubleStranded = False, outputDirectory = args.outputDirectory, clobber = self.mPileupAvatar.clobber)
        self.extractNeurosphere = programRunners.ExtractVariantsTumor(args.sampleName + ".neurosphere", self.mPileupNeurosphere.mPileupOut, minSupport = args.coverageRequirement, requireDoubleStranded = False, outputDirectory = args.outputDirectory, clobber = self.extractTumor.clobber)
        self.extractAvatar = programRunners.ExtractVariantsTumor(args.sampleName + ".avatar", self.mPileupAvatar.mPileupOut, minSupport = args.coverageRequirement, requireDoubleStranded = False, outputDirectory = args.outputDirectory, clobber = self.extractNeurosphere.clobber)
        #delete mPileups
        self.deleteNormalPileup = houseKeeping.Delete(self.mPileupNormal.mPileupOut)
        self.deleteTumorPileup = houseKeeping.Delete(self.mPileupTumor.mPileupOut)
        self.deleteNeurospherePileup = houseKeeping.Delete(self.mPileupNeurosphere.mPileupOut)
        self.deleteAvatarPileup = houseKeeping.Delete(self.mPileupAvatar.mPileupOut)
        #extract matched variant lines from the normal for comparisons
        self.extractNormalVsTumor = programRunners.ExtractVariantsNormal(args.sampleName + ".norm2tumor", self.mPileupNormal.mPileupOut, self.extractTumor.targetList, "normal2tumor", minSupport = args.coverageRequirement, clobber = self.extractAvatar.clobber, outputDirectory = args.outputDirectory)
        self.extractNormalVsNeurosphere = programRunners.ExtractVariantsNormal(args.sampleName + ".norm2neurosphere", self.mPileupNormal.mPileupOut, self.extractNeurosphere.targetList, "normal2neurosphere", minSupport = args.coverageRequirement, clobber = self.extractNormalVsTumor.clobber, outputDirectory = args.outputDirectory)
        self.extractNormalVsAvatar = programRunners.ExtractVariantsNormal(args.sampleName + ".norm2avatar", self.mPileupNormal.mPileupOut, self.extractAvatar.targetList, "normal2avatar", minSupport = args.coverageRequirement, clobber = self.extractNormalVsNeurosphere.clobber, outputDirectory = args.outputDirectory)
        #generate combines outputs
        self.combineNormalVsTumor = programRunners.CombineExtractedVariants(args.sampleName + ".tumorVcf", self.extractTumor.variantsOut, self.extractNormalVsTumor.variantsOut , "normalToTumor", self.extractNormalVsTumor.clobber, args.outputDirectory)
        self.combineNormalVsNeurosphere = programRunners.CombineExtractedVariants(args.sampleName + ".sphereVcf", self.extractNeurosphere.variantsOut, self.extractNormalVsNeurosphere.variantsOut , "normalToNeurosphere", self.combineNormalVsTumor.clobber, args.outputDirectory)
        self.combineNormalVsAvatar = programRunners.CombineExtractedVariants(args.sampleName + ".avatarVcf", self.extractAvatar.variantsOut, self.extractNormalVsAvatar.variantsOut , "normalToAvatar", self.combineNormalVsNeurosphere.clobber, args.outputDirectory)
        self.commandList = [self.alignNormal.pe1Command,
                            self.alignNormal.pe2Command,
                            self.alignNormal.samCommand,
                            self.alignTumor.pe1Command,
                            self.alignTumor.pe2Command,
                            self.alignTumor.samCommand,
                            self.alignNeurosphere.pe1Command,
                            self.alignNeurosphere.pe2Command,
                            self.alignNeurosphere.samCommand,
                            self.alignAvatar.pe1Command,
                            self.alignAvatar.pe2Command,
                            self.alignAvatar.samCommand,
                            self.makeBAMNormal.samToBamCommand,
                            self.makeBAMTumor.samToBamCommand,
                            self.makeBAMNeurosphere.samToBamCommand,
                            self.makeBAMAvatar.samToBamCommand,
                            self.deduplicateNormal.deduplicateCommand,
                            self.deduplicateTumor.deduplicateCommand,
                            self.deduplicateNeurosphere.deduplicateCommand,
                            self.deduplicateAvatar.deduplicateCommand,
                            self.mPileupNormal.mPileupCommand,
                            self.mPileupTumor.mPileupCommand,
                            self.mPileupNeurosphere.mPileupCommand,
                            self.mPileupAvatar.mPileupCommand,
                            self.extractTumor.extractCommand,
                            self.extractNeurosphere.extractCommand,
                            self.extractAvatar.extractCommand,
                            self.extractNormalVsTumor.extractCommand,
                            self.extractNormalVsNeurosphere.extractCommand,
                            self.extractNormalVsAvatar.extractCommand,
                            self.combineNormalVsTumor.combineCommand,
                            self.combineNormalVsNeurosphere.combineCommand,
                            self.combineNormalVsAvatar.combineCommand,
                            self.deleteNormalSai.deleteCommand,
                            self.deleteTumorSai.deleteCommand,
                            self.deleteNeurosphereSai.deleteCommand,
                            self.deleteAvatarSai.deleteCommand,
                            self.deleteNormalSam.deleteCommand,
                            self.deleteTumorSam.deleteCommand,
                            self.deleteNeurosphereSam.deleteCommand,
                            self.deleteAvatarSam.deleteCommand,
                            self.deleteFirstNormalBam.deleteCommand,
                            self.deleteFirstTumorBam.deleteCommand,
                            self.deleteFirstNeurosphereBam.deleteCommand,
                            self.deleteFirstAvatarBam.deleteCommand,
                            self.deleteNormalPileup.deleteCommand,
                            self.deleteTumorPileup.deleteCommand,
                            self.deleteNeurospherePileup.deleteCommand,
                            self.deleteAvatarPileup.deleteCommand]
        
def submitJobs(jobList):
    from runners import genericRunners
    tempDir = genericRunners.createTempDir(args.sampleName)
    writeCommandList(jobList.commandList, tempDir)
    alignNormalpe1 = genericRunners.HoffmanJob(1, [], args.sampleName + "na1", tempDir, args.emailAddress, "ba", cores = jobList.alignNormal.cores, mock = args.mock)
    alignNormalpe2 = genericRunners.HoffmanJob(2, [], args.sampleName + "na2", tempDir, args.emailAddress, args.emailConditions, cores = jobList.alignNormal.cores, mock = args.mock)
    alignNormalsam = genericRunners.HoffmanJob(3, [alignNormalpe1.jobID, alignNormalpe2.jobID], args.sampleName + "nsam", tempDir, args.emailAddress, args.emailConditions, mock = args.mock)
    alignTumorpe1 = genericRunners.HoffmanJob(4, [], args.sampleName + "ta1", tempDir, args.emailAddress, args.emailConditions, cores = jobList.alignTumor.cores, mock = args.mock)
    alignTumorpe2 = genericRunners.HoffmanJob(5, [], args.sampleName + "ta2", tempDir, args.emailAddress, args.emailConditions, cores = jobList.alignTumor.cores, mock = args.mock)
    alignTumorsam = genericRunners.HoffmanJob(6, [alignTumorpe1.jobID, alignTumorpe2.jobID], args.sampleName + "tsam", tempDir, args.emailAddress, args.emailConditions, mock = args.mock)
    alignNeurospherepe1 = genericRunners.HoffmanJob(7, [], args.sampleName + "sa1", tempDir, args.emailAddress, args.emailConditions, cores = jobList.alignNeurosphere.cores, mock = args.mock)
    alignNeurospherepe2 = genericRunners.HoffmanJob(8, [], args.sampleName + "sa2", tempDir, args.emailAddress, args.emailConditions, cores = jobList.alignNeurosphere.cores, mock = args.mock)
    alignNeurospheresam = genericRunners.HoffmanJob(9, [alignNeurospherepe1.jobID, alignNeurospherepe2.jobID], args.sampleName + "ssam", tempDir, args.emailAddress, args.emailConditions, mock = args.mock)
    alignAvatarpe1 = genericRunners.HoffmanJob(10, [], args.sampleName + "aa1", tempDir, args.emailAddress, args.emailConditions, cores = jobList.alignAvatar.cores, mock = args.mock)
    alignAvatarpe2 = genericRunners.HoffmanJob(11, [], args.sampleName + "aa2", tempDir, args.emailAddress, args.emailConditions, cores = jobList.alignAvatar.cores, mock = args.mock)
    alignAvatarsam = genericRunners.HoffmanJob(12, [alignAvatarpe1.jobID, alignAvatarpe2.jobID], args.sampleName + "asam", tempDir, args.emailAddress, args.emailConditions, mock = args.mock)
    makeBAMNormalsamToBam = genericRunners.HoffmanJob(13, [alignNormalsam.jobID], args.sampleName + "bn", tempDir, args.emailAddress, args.emailConditions, mock = args.mock)
    makeBAMTumorsamToBam = genericRunners.HoffmanJob(14, [alignTumorsam.jobID], args.sampleName + "bt", tempDir, args.emailAddress, args.emailConditions, mock = args.mock)
    makeBAMNeurospheresamToBam = genericRunners.HoffmanJob(15, [alignNeurospheresam.jobID], args.sampleName + "bs", tempDir, args.emailAddress, args.emailConditions, mock = args.mock)
    makeBAMAvatarsamToBam = genericRunners.HoffmanJob(16, [alignAvatarsam.jobID], args.sampleName + "ba", tempDir, args.emailAddress, args.emailConditions, mock = args.mock)
    deduplicateNormaldeduplicate = genericRunners.HoffmanJob(17, [makeBAMNormalsamToBam.jobID], args.sampleName + "dn", tempDir, args.emailAddress, args.emailConditions, mock = args.mock)
    deduplicateTumordeduplicate = genericRunners.HoffmanJob(18, [makeBAMTumorsamToBam.jobID], args.sampleName + "dt", tempDir, args.emailAddress, args.emailConditions, mock = args.mock)
    deduplicateNeurospherededuplicate = genericRunners.HoffmanJob(19, [makeBAMNeurospheresamToBam.jobID], args.sampleName + "ds", tempDir, args.emailAddress, args.emailConditions, mock = args.mock)
    deduplicateAvatardeduplicate = genericRunners.HoffmanJob(20, [makeBAMAvatarsamToBam.jobID], args.sampleName + "da", tempDir, args.emailAddress, args.emailConditions, mock = args.mock)
    mPileupNormalmPileup = genericRunners.HoffmanJob(21, [deduplicateNormaldeduplicate.jobID], args.sampleName + "mn", tempDir, args.emailAddress, args.emailConditions, mock = args.mock)
    mPileupTumormPileup = genericRunners.HoffmanJob(22, [deduplicateTumordeduplicate.jobID], args.sampleName + "mt", tempDir, args.emailAddress, args.emailConditions, mock = args.mock)
    mPileupNeurospheremPileup = genericRunners.HoffmanJob(23, [deduplicateNeurospherededuplicate.jobID], args.sampleName + "ms", tempDir, args.emailAddress, args.emailConditions, mock = args.mock)
    mPileupAvatarmPileup = genericRunners.HoffmanJob(24, [deduplicateAvatardeduplicate.jobID], args.sampleName + "ma", tempDir, args.emailAddress, args.emailConditions, mock = args.mock)
    extractTumorextract = genericRunners.HoffmanJob(25, [mPileupTumormPileup.jobID], args.sampleName + "et", tempDir, args.emailAddress, args.emailConditions, mock = args.mock)
    extractNeurosphereextract = genericRunners.HoffmanJob(26, [mPileupNeurospheremPileup.jobID], args.sampleName + "en", tempDir, args.emailAddress, args.emailConditions, mock = args.mock)
    extractAvatarextract = genericRunners.HoffmanJob(27, [mPileupAvatarmPileup.jobID], args.sampleName + "ea", tempDir, args.emailAddress, args.emailConditions, mock = args.mock)
    extractNormalVsTumorextract = genericRunners.HoffmanJob(28, [mPileupNormalmPileup.jobID, extractTumorextract.jobID], args.sampleName + "nvt", tempDir, args.emailAddress, args.emailConditions, mock = args.mock)
    extractNormalVsNeurosphereextract = genericRunners.HoffmanJob(29, [mPileupNormalmPileup.jobID, extractNeurosphereextract.jobID], args.sampleName + "nvs", tempDir, args.emailAddress, args.emailConditions, mock = args.mock)
    extractNormalVsAvatarextract = genericRunners.HoffmanJob(30, [mPileupNormalmPileup.jobID, extractAvatarextract.jobID], args.sampleName + "nva", tempDir, args.emailAddress, args.emailConditions, mock = args.mock)
    combineTumorVsNormal = genericRunners.HoffmanJob(31, [extractTumorextract.jobID, extractNormalVsTumorextract.jobID], args.sampleName + "tvcf", tempDir, args.emailAddress, args.emailConditions, mock = args.mock)
    combineNeurosphereVsNormal = genericRunners.HoffmanJob(32, [extractNeurosphereextract.jobID, extractNormalVsNeurosphereextract.jobID], args.sampleName + "nvcf", tempDir, args.emailAddress, args.emailConditions, mock = args.mock)
    combineAvatarVsNormal = genericRunners.HoffmanJob(33, [extractAvatarextract.jobID, extractNormalVsAvatarextract.jobID], args.sampleName + "avcf", tempDir, args.emailAddress, "ae", mock = args.mock)
    deleteNormalSai = genericRunners.HoffmanJob(34, [alignNormalsam.jobID], args.sampleName + "delnsai", tempDir, args.emailAddress, args.emailConditions, mock = args.mock)
    deleteTumorSai = genericRunners.HoffmanJob(35, [alignTumorsam.jobID], args.sampleName + "deltsai", tempDir, args.emailAddress, args.emailConditions, mock = args.mock)
    deleteNeurosphereSai = genericRunners.HoffmanJob(36, [alignNeurospheresam.jobID], args.sampleName + "delssai", tempDir, args.emailAddress, args.emailConditions, mock = args.mock)
    deleteAvatarSai = genericRunners.HoffmanJob(37, [alignAvatarsam.jobID], args.sampleName + "delasai", tempDir, args.emailAddress, args.emailConditions, mock = args.mock)
    deleteNormalSam = genericRunners.HoffmanJob(38, [makeBAMNormalsamToBam.jobID], args.sampleName + "delnsam", tempDir, args.emailAddress, args.emailConditions, mock = args.mock)
    deleteTumorSam = genericRunners.HoffmanJob(39, [makeBAMTumorsamToBam.jobID], args.sampleName + "deltsam", tempDir, args.emailAddress, args.emailConditions, mock = args.mock)
    deleteNeurosphereSam = genericRunners.HoffmanJob(40, [makeBAMNeurospheresamToBam.jobID], args.sampleName + "delssam", tempDir, args.emailAddress, args.emailConditions, mock = args.mock)
    deleteAvatarSam = genericRunners.HoffmanJob(41, [makeBAMAvatarsamToBam.jobID], args.sampleName + "delasam", tempDir, args.emailAddress, args.emailConditions, mock = args.mock)
    deleteFirstNormalBam = genericRunners.HoffmanJob(42, [deduplicateNormaldeduplicate.jobID], args.sampleName + "delnbam", tempDir, args.emailAddress, args.emailConditions, mock = args.mock)
    deleteFirstTumorBam = genericRunners.HoffmanJob(43, [deduplicateTumordeduplicate.jobID], args.sampleName + "deltbam", tempDir, args.emailAddress, args.emailConditions, mock = args.mock)
    deleteFirstNeurosphereBam = genericRunners.HoffmanJob(44, [deduplicateNeurospherededuplicate.jobID], args.sampleName + "delsbam", tempDir, args.emailAddress, args.emailConditions, mock = args.mock)
    deleteFirstAvatarBam = genericRunners.HoffmanJob(45, [deduplicateAvatardeduplicate.jobID], args.sampleName + "delbbam", tempDir, args.emailAddress, args.emailConditions, mock = args.mock)
    deleteNormalPileup = genericRunners.HoffmanJob(46, [extractNormalVsAvatarextract.jobID, extractNormalVsNeurosphereextract.jobID, extractNormalVsTumorextract.jobID], args.sampleName + "delnmp", tempDir, args.emailAddress, args.emailConditions, mock = args.mock)
    deleteTumorPileup = genericRunners.HoffmanJob(47, [extractTumorextract.jobID, extractNormalVsTumorextract.jobID], args.sampleName + "deltmp", tempDir, args.emailAddress, args.emailConditions, mock = args.mock)
    deleteNeurospherePileup = genericRunners.HoffmanJob(48, [extractNeurosphereextract.jobID, extractNormalVsNeurosphereextract.jobID], args.sampleName + "delsmp", tempDir, args.emailAddress, args.emailConditions, mock = args.mock)
    deleteAvatarPileup = genericRunners.HoffmanJob(49, [extractAvatarextract.jobID, extractNormalVsAvatarextract.jobID], args.sampleName + "delamp", tempDir, args.emailAddress, args.emailConditions, mock = args.mock)
    
    
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