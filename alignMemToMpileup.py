#!/usr/bin/env python3

class CheckArgs(object):
    
    def __init__(self):
        import argparse
        import os
        parser = argparse.ArgumentParser()
        parser.add_argument("-m", "--mock", help = "Mock submission (don't actually submit to scheduler)", action = 'store_true')
        parser.add_argument("-9", "--clobber", help = "Overwrite existing files without asking.", action = 'store_true')
        parser.add_argument("-1", "--pe1", help = "Paired-end 1 fastq file.")
        parser.add_argument("-2", "--pe2", help = "Paired-end 2 fastq file.", default = "")
        parser.add_argument("-n", "--sampleName", help = "The unique sample identifier. Seriously, make it unique.  Or not... be lazy and reap the consequences.")
        parser.add_argument("-r", "--referenceGenome", help = "The reference genome fasta")
        parser.add_argument("-M", "--emailAddress", help = "Email address for notifications", default = "")
        parser.add_argument("--emailConditions", help = "Conditions to email for (b)egin, (e)nd, (a)bort.", default = "a")
        rawArgs = parser.parse_args()
        self.mock = rawArgs.mock
        clobber = rawArgs.clobber
        if not clobber:
            self.clobber = None  #using none as a default type here.  We will force the user to set this if it becomes an issue.  Hopefully it won't because they are using unique names.
        else:
            self.clobber = True
        pe1 = rawArgs.pe1
        if not pe1:
            raise RuntimeError("No file specified for paired-end 1.  Nothing to analyze here.")
        if not os.path.isfile(pe1):
            raise FileNotFoundError("Paired-end 1 file was not found. Expected file: %s" %(pe1))
        self.pe1 = pe1
        if not rawArgs.pe2:
            self.pe2 = ""
        else:
            pe2 = rawArgs.pe2
            if not os.path.isfile(pe2):
                raise FileNotFoundError("Paired-end 2 file was not found. Expected file: %s" %(pe2))
            self.pe2 = pe2
        if self.pe1 == self.pe2:
            raise RuntimeError("Both paired end files are the same file. This is most likely an error.")
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
            
class JobList():
    def __init__(self):
        from runners import programRunners
        self.alignment = programRunners.BWAmem(args.sampleName, args.referenceGenome, args.pe1, args.pe2, clobber = args.clobber)
        self.makeBAM = programRunners.SAMtoBAM(args.sampleName, self.alignment.samOut, clobber = self.alignment.clobber)  #clobber can be changed during the previous object initialization.  Using that to keep current on changes.
        self.deduplicate = programRunners.Deduplicate(args.sampleName, self.makeBAM.bamOut, clobber = self.makeBAM.clobber)
        self.mPileup = programRunners.MPileup(args.sampleName, self.deduplicate.bamOut, args.referenceGenome, clobber = self.deduplicate.clobber)
        self.pileupToVCF = programRunners.PileupToVCF(args.sampleName, self.mPileup.mPileupOut, clobber = self.mPileup.clobber)
        self.commandList = [self.alignment.bwaCommand,
                            self.makeBAM.samToBamCommand,
                            self.deduplicate.deduplicateCommand,
                            self.mPileup.mPileupCommand,
                            self.pileupToVCF.pileupToVCFCommand]
        
def submitJobs(jobList):
    from runners import genericRunners
    tempDir = genericRunners.createTempDir(args.sampleName)
    writeCommandList(jobList.commandList, tempDir)
    align = genericRunners.HoffmanJob([], 1, args.sampleName + "align", tempDir, args.emailAddress, "ba", cores = jobList.alignment.threads, mock = args.mock)
    makeBAM = genericRunners.HoffmanJob([makeSAM.jobID], 2, args.sampleName + "bam", tempDir, args.emailAddress, args.emailConditions, mock = args.mock)
    deduplicate = genericRunners.HoffmanJob([makeBAM.jobID], 3, args.sampleName + "dup", tempDir, args.emailAddress, args.emailConditions, mock = args.mock)
    runMPileup = genericRunners.HoffmanJob([makeBAM.jobID], 4, args.sampleName + "pileup", tempDir, args.emailAddress, args.emailConditions, mock = args.mock)
    runPileupToVCF = genericRunners.HoffmanJob([runMPileup.jobID], 5, args.sampleName + "vcf", tempDir, args.emailAddress, "ea", mock = args.mock)
    
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
    jobList = JobList()
    submitJobs(jobList)
    print("Completed")
    
main()