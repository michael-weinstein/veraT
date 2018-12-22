#!/usr/bin/env python3

import os

class PipelineStart(object):
    
    def __init__(self, sampleName, useDir = False, tempDir = False):
        if tempDir:
            import os
            if not os.path.isdir(tempDir):
                os.mkdir(tempDir)
            else:
                if os.path.isfile(tempDir + os.sep + "jobs.pkl"):
                    raise RuntimeError("A jobs.pkl already exists in the specified temporary directory.  Please clear the directory before overwriting with the new jobs.")
        elif useDir:
            import runnerSupport
            import os
            if not os.path.isdir(useDir):
                os.mkdir(useDir)
            tempDir = runnerSupport.createTempDir(sampleName, useDir)
        else:
            import runnerSupport
            tempDir = runnerSupport.createTempDir(sampleName)
        self.tempDir = tempDir
        self.jobs = runnerSupport.JobList(self.tempDir)
        
    def getJobList(self):
        return self.jobs
        
    def end(self):
        self.jobs.close()
        
class AlignMemAddReadGroups(object):
    
    def __init__(self, jobList, sampleName, refGenomeFasta, pe1 = False, pe2 = False, bwaCores = "calculateFromSize", readGroupLibrary = "defaultLibrary", readGroupID = 1, readGroupSampleName = False, readGroupPlatform = "Illumina", barcode = "defaultBarcode", emailAddress = False, clobber = None, outputDir = "", lastPe1Job = False, lastPe2Job = False, mock = False):
        import programRunners
        import houseKeeping
        import picardRunners
        import runnerSupport
        self.emailAddress = emailAddress
        self.tempDir = jobList.tempDir
        self.sampleName = sampleName
        self.library = str(readGroupLibrary)
        if self.library == "defaultLibrary":
            self.library = ""
        if bwaCores == "calculateFromSize":
            self.bwaCores = runnerSupport.calculateCoresFromFileSize([pe1, pe2], 5)
        elif type(bwaCores) == int:
            self.bwaCores = bwaCores
        elif type(bwaCores) == float:
            self.bwaCores = max([int(bwaCores), 1])
        else:
            raise RuntimeError("BWA cores argument was not a number and not a direction to calculation from size.  Passed value: %s" %(bwaCores))
        if not lastPe1Job and not pe1:
            raise RuntimeError("No paired-end 1 file set (nothing to align or pe2 set without a pe1)")
        if not pe1:
            pe1 = lastPe1Job.data
        if lastPe2Job and not pe2:
            pe2 = lastPe2Job.data
        if not lastPe1Job:
            lastPe1Job = runnerSupport.EmptyWorkflowReturn()
        if not lastPe2Job:
            lastPe2Job = runnerSupport.EmptyWorkflowReturn()
        self.lastPe1Job = lastPe1Job
        self.lastPe2Job = lastPe2Job
        clobberSet = set([clobber] + [item.clobber for item in [lastPe1Job, lastPe2Job]])
        if True in clobberSet:
            clobber = True
        elif False  in clobberSet:
            clobber = False
        elif None in clobberSet:
            clobber = None
        align = programRunners.BWAmem(sampleName, refGenomeFasta, pe1, pe2, cores = self.bwaCores, clobber = None, outputDirectory = outputDir)
        self.bwaCommandIndex = jobList.addJob(align.bwaCommand)
        self.samFile = align.samOut
        makeBam = programRunners.ViewSAMtoBAM(sampleName, self.samFile, refGenomeFasta, clobber = align.clobber, outputDirectory = outputDir)
        self.viewSamCommandIndex = jobList.addJob(makeBam.viewCommand)
        self.bamFile = makeBam.bamOut
        deleteSam = houseKeeping.Delete(self.samFile)
        self.deleteSamCommandIndex = jobList.addJob(deleteSam.deleteCommand)
        addReadGroups = picardRunners.AddReadGroups(sampleName, self.bamFile, rgid = readGroupID, rglb = readGroupLibrary, rgpl = readGroupPlatform, rgsm = readGroupSampleName, rgpu = barcode, clobber = makeBam.clobber, outputDirectory = outputDir)
        self.addReadGroupsCommandIndex = jobList.addJob(addReadGroups.addReadGroupsCommand)
        deleteOriginalBAM = houseKeeping.Delete(makeBam.bamOut)
        self.deleteOriginalBAMCommandIndex = jobList.addJob(deleteOriginalBAM.deleteCommand)
        addReadGroupsJobID = self.submitCommands(mock)
        self.returnData = runnerSupport.WorkflowReturn(addReadGroups.bamOut, addReadGroupsJobID, addReadGroups.clobber)
        
    def submitCommands(self, mock):
        import genericRunners
        tempDir = self.tempDir
        if self.bwaCommandIndex == 1:
            emailStatus = "ba"
        else:
            emailStatus = "a"
        align = genericRunners.HoffmanJob([self.lastPe1Job.jobID, self.lastPe2Job.jobID], self.bwaCommandIndex, self.sampleName + self.library + "align", tempDir, self.emailAddress, emailStatus, cores = self.bwaCores, mock = mock)
        encodeBAM = genericRunners.HoffmanJob([align.jobID], self.viewSamCommandIndex, self.sampleName + self.library + "BAM", tempDir, self.emailAddress, "a", cores = 1, mock = mock)
        deleteSAM = genericRunners.HoffmanJob([encodeBAM.jobID], self.deleteSamCommandIndex, self.sampleName + self.library + "delSAM", tempDir, self.emailAddress, "a", cores = 1, mock = mock)
        addReadGroups = genericRunners.HoffmanJob([encodeBAM.jobID], self.addReadGroupsCommandIndex, self.sampleName + self.library + "AddRG", tempDir, self.emailAddress, "a", cores = 1, mock = mock)
        deleteOriginalBAM = genericRunners.HoffmanJob([addReadGroups.jobID], self.deleteOriginalBAMCommandIndex, self.sampleName + self.library + "delOrigBAM", tempDir, self.emailAddress, "a", cores = 1, mock = mock)
        return addReadGroups.jobID
    
class AlignAlnAddReadGroups(object):
    
    def __init__(self, jobList, sampleName, refGenomeFasta, pe1 = False, pe2 = False, bwaCores = "calculateFromSize", readGroupLibrary = "defaultLibrary", readGroupID = 1, readGroupSampleName = False, readGroupPlatform = "Illumina", barcode = "defaultBarcode", emailAddress = False, flagFilterOut = False, clobber = None, outputDir = "", lastPe1Job = False, lastPe2Job = False, mock = False):
        import programRunners
        import houseKeeping
        import picardRunners
        import runnerSupport
        self.emailAddress = emailAddress
        self.tempDir = jobList.tempDir
        self.sampleName = sampleName
        self.library = str(readGroupLibrary)
        if self.library == "defaultLibrary":
            self.library = ""
        if bwaCores == "calculateFromSize":
            self.bwaCores = runnerSupport.calculateCoresFromFileSize([pe1, pe2], 5)
        elif type(bwaCores) == int:
            self.bwaCores = bwaCores
        elif type(bwaCores) == float:
            self.bwaCores = max([int(bwaCores), 1])
        else:
            raise RuntimeError("BWA cores argument was not a number and not a direction to calculation from size.  Passed value: %s" %(bwaCores))
        if not lastPe1Job and not pe1:
            raise RuntimeError("No paired-end 1 file set (nothing to align or pe2 set without a pe1)")
        if not pe1:
            pe1 = lastPe1Job.data
        if lastPe2Job and not pe2:
            pe2 = lastPe2Job.data
        if not lastPe1Job:
            lastPe1Job = runnerSupport.EmptyWorkflowReturn()
        if not lastPe2Job:
            lastPe2Job = runnerSupport.EmptyWorkflowReturn()
        self.lastPe1Job = lastPe1Job
        self.lastPe2Job = lastPe2Job
        clobberSet = set([clobber] + [item.clobber for item in [lastPe1Job, lastPe2Job]])
        if True in clobberSet:
            clobber = True
        elif False  in clobberSet:
            clobber = False
        elif None in clobberSet:
            clobber = None
        align = programRunners.BWAlign(sampleName, refGenomeFasta, pe1, pe2, cores = self.bwaCores, clobber = clobber, outputDirectory = outputDir)
        self.alignPe1CommandIndex = jobList.addJob(align.pe1Command)
        if align.pe2Command:
            self.alignPe2CommandIndex = jobList.addJob(align.pe2Command)
        else:
            self.alignPe2CommandIndex = False
        self.makeSamCommandIndex = jobList.addJob(align.samCommand)
        self.samFile = align.samOut
        makeBam = programRunners.ViewSAMtoBAM(sampleName, self.samFile, refGenomeFasta, flagFilterOut = flagFilterOut, clobber = align.clobber, outputDirectory = outputDir)
        deleteSai1 = houseKeeping.Delete(align.pe1Out)
        self.deleteSai1CommandIndex = jobList.addJob(deleteSai1.deleteCommand)
        if align.pe2Out:
            deleteSai2 = houseKeeping.Delete(align.pe2Out)
            self.deleteSai2CommandIndex = jobList.addJob(deleteSai2.deleteCommand)
        else:
            self.deleteSai2CommandIndex = False
        self.viewSamCommandIndex = jobList.addJob(makeBam.viewCommand)
        self.bamFile = makeBam.bamOut
        deleteSam = houseKeeping.Delete(self.samFile)
        self.deleteSamCommandIndex = jobList.addJob(deleteSam.deleteCommand)
        addReadGroups = picardRunners.AddReadGroups(sampleName, self.bamFile, rgid = readGroupID, rglb = readGroupLibrary, rgpl = readGroupPlatform, rgsm = readGroupSampleName, rgpu = barcode, clobber = makeBam.clobber, outputDirectory = outputDir)
        self.addReadGroupsCommandIndex = jobList.addJob(addReadGroups.addReadGroupsCommand)
        deleteOriginalBAM = houseKeeping.Delete(makeBam.bamOut)
        self.deleteOriginalBAMCommandIndex = jobList.addJob(deleteOriginalBAM.deleteCommand)
        addReadGroupsJobID = self.submitCommands(mock)
        self.returnData = runnerSupport.WorkflowReturn(addReadGroups.bamOut, addReadGroupsJobID, addReadGroups.clobber)
        
    def submitCommands(self, mock):
        import genericRunners
        import runnerSupport
        tempDir = self.tempDir
        if self.alignPe1CommandIndex == 1:
            emailStatus = "ba"
        else:
            emailStatus = "a"
        align1 = genericRunners.HoffmanJob([self.lastPe1Job.jobID], self.alignPe1CommandIndex, self.sampleName + self.library + "align1", tempDir, self.emailAddress, emailStatus, cores = self.bwaCores, mock = mock)
        if self.alignPe2CommandIndex:
            align2 = genericRunners.HoffmanJob([self.lastPe2Job.jobID], self.alignPe2CommandIndex, self.sampleName + self.library + "align2", tempDir, self.emailAddress, emailStatus, cores = self.bwaCores, mock = mock)
        else:
            align2 = runnerSupport.EmptyWorkflowReturn()
        makeSam = genericRunners.HoffmanJob([align1.jobID, align2.jobID], self.makeSamCommandIndex, self.sampleName + self.library + "SAM", tempDir, self.emailAddress, emailStatus, cores = 1, mock = mock)
        deleteSai1 = genericRunners.HoffmanJob([makeSam.jobID], self.deleteSai1CommandIndex, self.sampleName + self.library + "delSai1", tempDir, self.emailAddress, emailStatus, cores = 1, mock = mock)
        if self.deleteSai2CommandIndex:
            deleteSai2 = genericRunners.HoffmanJob([makeSam.jobID], self.deleteSai2CommandIndex, self.sampleName + self.library + "delSai2", tempDir, self.emailAddress, emailStatus, cores = 1, mock = mock)
        encodeBAM = genericRunners.HoffmanJob([makeSam.jobID], self.viewSamCommandIndex, self.sampleName + self.library + "BAM", tempDir, self.emailAddress, "a", cores = 1, mock = mock)
        deleteSAM = genericRunners.HoffmanJob([encodeBAM.jobID], self.deleteSamCommandIndex, self.sampleName + self.library + "delSAM", tempDir, self.emailAddress, "a", cores = 1, mock = mock)
        addReadGroups = genericRunners.HoffmanJob([encodeBAM.jobID], self.addReadGroupsCommandIndex, self.sampleName + self.library + "AddRG", tempDir, self.emailAddress, "a", cores = 1, mock = mock)
        deleteOriginalBAM = genericRunners.HoffmanJob([addReadGroups.jobID], self.deleteOriginalBAMCommandIndex, self.sampleName + self.library + "delOrigBAM", tempDir, self.emailAddress, "a", cores = 1, mock = mock)
        return addReadGroups.jobID
    
class AlignHisat2AddReadGroups(object):
    
    def __init__(self, jobList, sampleName, refGenomeFasta, refGenomePrefix, pe1 = False, pe2 = False, hisatCores = "calculateFromSize", readGroupLibrary = "defaultLibrary", readGroupID = 1, readGroupSampleName = False, readGroupPlatform = "Illumina", barcode = "defaultBarcode", emailAddress = False, clobber = None, outputDir = "", lastPe1Job = False, lastPe2Job = False, mock = False):
        import programRunners
        import houseKeeping
        import picardRunners
        import runnerSupport
        self.emailAddress = emailAddress
        self.tempDir = jobList.tempDir
        self.sampleName = sampleName
        self.library = str(readGroupLibrary)
        if self.library == "defaultLibrary":
            self.library = ""
        if hisatCores == "calculateFromSize":
            self.hisatCores = runnerSupport.calculateCoresFromFileSize([pe1, pe2], 5)
        elif type(hisatCores) == int:
            self.hisatCores = hisatCores
        elif type(hisatCores) == float:
            self.hisatCores = max([int(hisatCores), 1])
        else:
            raise RuntimeError("BWA cores argument was not a number and not a direction to calculation from size.  Passed value: %s" %(hisatCores))
        if not lastPe1Job and not pe1:
            raise RuntimeError("No paired-end 1 file set (nothing to align or pe2 set without a pe1)")
        if not pe1:
            pe1 = lastPe1Job.data
        if lastPe2Job and not pe2:
            pe2 = lastPe2Job.data
        if not lastPe1Job:
            lastPe1Job = runnerSupport.EmptyWorkflowReturn()
        if not lastPe2Job:
            lastPe2Job = runnerSupport.EmptyWorkflowReturn()
        self.lastPe1Job = lastPe1Job
        self.lastPe2Job = lastPe2Job
        clobberSet = set([clobber] + [item.clobber for item in [lastPe1Job, lastPe2Job]])
        if True in clobberSet:
            clobber = True
        elif False  in clobberSet:
            clobber = False
        elif None in clobberSet:
            clobber = None
        align = programRunners.Hisat2Align(sampleName, refGenomePrefix, pe1, pe2, cores = self.hisatCores, clobber = None, outputDirectory = outputDir)
        self.hisat2CommandIndex = jobList.addJob(align.hisat2Command)
        self.samFile = align.samOut
        makeBam = programRunners.ViewSAMtoBAM(sampleName, self.samFile, refGenomeFasta, clobber = align.clobber, outputDirectory = outputDir)
        self.viewSamCommandIndex = jobList.addJob(makeBam.viewCommand)
        self.bamFile = makeBam.bamOut
        deleteSam = houseKeeping.Delete(self.samFile)
        self.deleteSamCommandIndex = jobList.addJob(deleteSam.deleteCommand)
        addReadGroups = picardRunners.AddReadGroups(sampleName, self.bamFile, rgid = readGroupID, rglb = readGroupLibrary, rgpl = readGroupPlatform, rgsm = readGroupSampleName, rgpu = barcode, clobber = makeBam.clobber, outputDirectory = outputDir)
        self.addReadGroupsCommandIndex = jobList.addJob(addReadGroups.addReadGroupsCommand)
        deleteOriginalBAM = houseKeeping.Delete(makeBam.bamOut)
        self.deleteOriginalBAMCommandIndex = jobList.addJob(deleteOriginalBAM.deleteCommand)
        addReadGroupsJobID = self.submitCommands(mock)
        self.returnData = runnerSupport.WorkflowReturn(addReadGroups.bamOut, addReadGroupsJobID, addReadGroups.clobber)
        
    def submitCommands(self, mock):
        import genericRunners
        tempDir = self.tempDir
        if self.hisat2CommandIndex == 1:
            emailStatus = "ba"
        else:
            emailStatus = "a"
        align = genericRunners.HoffmanJob([], self.hisat2CommandIndex, self.sampleName + self.library + "align", tempDir, self.emailAddress, emailStatus, cores = self.hisatCores, mock = mock)
        encodeBAM = genericRunners.HoffmanJob([align.jobID], self.viewSamCommandIndex, self.sampleName + self.library + "BAM", tempDir, self.emailAddress, "a", cores = 1, mock = mock)
        deleteSAM = genericRunners.HoffmanJob([encodeBAM.jobID], self.deleteSamCommandIndex, self.sampleName + self.library + "delSAM", tempDir, self.emailAddress, "a", cores = 1, mock = mock)
        addReadGroups = genericRunners.HoffmanJob([encodeBAM.jobID], self.addReadGroupsCommandIndex, self.sampleName + self.library + "AddRG", tempDir, self.emailAddress, "a", cores = 1, mock = mock)
        deleteOriginalBAM = genericRunners.HoffmanJob([addReadGroups.jobID], self.deleteOriginalBAMCommandIndex, self.sampleName + self.library + "delOrigBAM", tempDir, self.emailAddress, "a", cores = 1, mock = mock)
        return addReadGroups.jobID
    
class MergeBAMFilesAndDeduplicate(object):
    
    def __init__(self, jobList, sampleName, bamFiles = False, sort_order = "coordinate", validation_stringency = "LENIENT", create_index = True, emailAddress = False, clobber = None, outputDir = "", lastJob = False, mock = False):
        import houseKeeping
        import picardRunners
        import runnerSupport
        if not bamFiles and not lastJob:
            raise RuntimeError("No bams for processing specified in last job data or arguments.")
        if not lastJob:
            lastJob = runnerSupport.EmptyWorkflowReturn()
        else:
            clobber = lastJob.clobber
        self.lastJob = lastJob
        if not bamFiles:
            bamFiles = lastJob.data
        if type(bamFiles) == str:
            bamFiles = [bamFiles]
        self.emailAddress = emailAddress
        self.tempDir = jobList.tempDir
        self.sampleName = sampleName
        self.mergeCommandIndex = False  #default value to be changed if a merge is required
        if len(bamFiles) > 1:  #if we have more than one bam file come in, we need to merge them first, otherwise we can skip the merge process
            merge = picardRunners.MergeSAMFiles(sampleName, bamFiles, sort_order, validation_stringency, create_index, clobber, outputDir)
            self.mergeCommandIndex = jobList.addJob(merge.mergeSAMCommand)
            self.deleteCommands = []
            for bamFile in bamFiles:
                deletePartialBAM = houseKeeping.Delete(bamFile)
                self.deleteCommands.append(deletePartialBAM.deleteCommand)
            self.deletePartialCommandIndex = jobList.addJob(self.deleteCommands)  #going to have just one node do all the deletions for this job
        if len(bamFiles) > 1:  #if we have multiple bam files and had to merge them
            deduplicate = picardRunners.Deduplicate(sampleName, merge.bamOut, clobber = merge.clobber, outputDirectory = outputDir)
        else:  #if we had one bam file come in, we can skip merging and go directly to deduplication
            deduplicate = picardRunners.Deduplicate(sampleName, bamFiles[0], revertToBaseName = True, clobber = clobber, outputDirectory = outputDir)
        self.deduplicateCommandIndex = jobList.addJob(deduplicate.deduplicateCommand)
        if len(bamFiles) > 1:
            deleteOriginalBam = houseKeeping.Delete(merge.bamOut)
        else:
            deleteOriginalBam = houseKeeping.Delete(bamFiles[0])
        self.deleteOriginalBAMCommandIndex = jobList.addJob(deleteOriginalBam.deleteCommand)
        deduplicateJobID = self.submitCommands(mock)
        self.returnData = runnerSupport.WorkflowReturn(deduplicate.bamOut, deduplicateJobID, deduplicate.clobber)
        
    def submitCommands(self, mock):
        import genericRunners
        tempDir = self.tempDir
        if self.mergeCommandIndex:
            merge = genericRunners.HoffmanJob(self.lastJob.jobID, self.mergeCommandIndex, self.sampleName + "Merge", tempDir, self.emailAddress, "a", 1, mock = mock)
            deletePartialBAM = genericRunners.HoffmanJob(merge.jobID, self.deletePartialCommandIndex, self.sampleName + "DelPart", tempDir, self.emailAddress, "a", cores = 1, mock = mock)
            dedupeHold = merge.jobID
        else:
            dedupeHold = self.lastJob.jobID
        deduplicate = genericRunners.HoffmanJob(dedupeHold, self.deduplicateCommandIndex, self.sampleName + "Dedup", tempDir, self.emailAddress, "a", cores = 1, mock = mock)
        deleteOriginalBAM = genericRunners.HoffmanJob(deduplicate.jobID, self.deleteOriginalBAMCommandIndex, self.sampleName + "DelOrig", tempDir, self.emailAddress, "a", cores = 1, mock = mock)
        return deduplicate.jobID
    
class GATKProcessBAM(object):
    
    def __init__(self, jobList, sampleName, refGenomeFasta, bamFile = False, dbSNP = False, intervals = False, emailAddress = False, clobber = None, outputDir = "", lastJob = False, mock = False):
        import houseKeeping
        import gatkRunners
        import runnerSupport
        import picardRunners
        if not bamFile and not lastJob:
            raise RuntimeError("No bam for processing specified in last job data or arguments.")
        if not lastJob:
            lastJob = runnerSupport.EmptyWorkflowReturn()
        else:
            clobber = lastJob.clobber
        self.lastJob = lastJob
        if not bamFile:
            bamFile = lastJob.data
        self.emailAddress = emailAddress
        self.tempDir = jobList.tempDir
        self.sampleName = sampleName
        indelRealignment = gatkRunners.IndelRealignment(sampleName, bamFile, refGenomeFasta, dbSNP, intervals, clobber = clobber, outputDirectory = outputDir)
        self.findTargetsCommandIndex = jobList.addJob(indelRealignment.targetCreationCommand)
        self.realignerCommandIndex = jobList.addJob(indelRealignment.realignerCommand)
        bqsr = gatkRunners.BQSR(sampleName, indelRealignment.bamOut, refGenomeFasta, dbSNP, intervals, clobber = indelRealignment.clobber, outputDirectory = outputDir)
        self.bqsrFirstPassCommandIndex = jobList.addJob(bqsr.firstPassCommand)
        self.bqsrExecuteCommandIndex = jobList.addJob(bqsr.printReadsCommand)
        self.bqsrAnalysisCommandIndex = jobList.addJob(bqsr.secondPassCommand)#, bqsr.plotCommand])
        indexRecalBAM = picardRunners.BuildBAMIndex(bqsr.bamOut)
        self.indexRecalBAMCommandIndex = jobList.addJob(indexRecalBAM.indexCommand)
        deleteRealignedBAM = houseKeeping.Delete(indelRealignment.bamOut)
        self.deleteRealignedBAMCommandIndex = jobList.addJob(deleteRealignedBAM.deleteCommand)
        depthOfCoverage = gatkRunners.DepthOfCoverage(sampleName, bqsr.bamOut, refGenomeFasta, intervals, clobber = bqsr.clobber, outputDirectory = outputDir)
        self.depthOfCoverageCommandIndex = jobList.addJob(depthOfCoverage.depthOfCoverageCommand)
        indexRecalBAMJobID = self.submitCommands(mock)
        self.returnData = runnerSupport.WorkflowReturn(bqsr.bamOut, indexRecalBAMJobID, depthOfCoverage.clobber)
        
    def submitCommands(self, mock):
        import genericRunners
        tempDir = self.tempDir
        findTargets = genericRunners.HoffmanJob(self.lastJob.jobID, self.findTargetsCommandIndex, self.sampleName + "indel1", tempDir, self.emailAddress, "a", 1, memory = 16, mock = mock)
        realigner = genericRunners.HoffmanJob(findTargets.jobID, self.realignerCommandIndex, self.sampleName + "indel2", tempDir, self.emailAddress, "a", 1, memory = 16, mock = mock)
        bqsrFirstPass = genericRunners.HoffmanJob(realigner.jobID, self.bqsrFirstPassCommandIndex, self.sampleName + "bqsr1", tempDir, self.emailAddress, "a", 1, memory = 16, mock = mock)
        printReads = genericRunners.HoffmanJob(bqsrFirstPass.jobID, self.bqsrExecuteCommandIndex, self.sampleName + "bqsr2", tempDir, self.emailAddress, "a", 1, memory = 16, mock = mock)
        indexRecalBAM = genericRunners.HoffmanJob(printReads.jobID, self.indexRecalBAMCommandIndex, self.sampleName + "idxRecal", tempDir, self.emailAddress, "a", 1, memory = 16, mock = mock)
        bqsrAnalysis = genericRunners.HoffmanJob(bqsrFirstPass.jobID, self.bqsrAnalysisCommandIndex, self.sampleName + "bqsrData", tempDir, self.emailAddress, "a", 1, memory = 16, mock = mock)
        deleteRealignedBAM = genericRunners.HoffmanJob([printReads.jobID, bqsrAnalysis.jobID], self.deleteRealignedBAMCommandIndex, self.sampleName + "DelRealign", tempDir, self.emailAddress, "a", 1, mock = mock)
        depthOfCoverage = genericRunners.HoffmanJob([printReads.jobID, indexRecalBAM.jobID], self.depthOfCoverageCommandIndex, self.sampleName + "depth", tempDir, self.emailAddress, "a", 1, memory = 16, mock = mock)
        return indexRecalBAM.jobID
    
class SplitAndTrimRNAData(object):
    
    def __init__(self, jobList, sampleName, refGenomeFasta, bamFile = False, emailAddress = False, clobber = None, outputDir = "", lastJob = False, mock = False):
        import houseKeeping
        import gatkRunners
        import runnerSupport
        if not bamFile and not lastJob:
            raise RuntimeError("No bam for processing specified in last job data or arguments.")
        if not lastJob:
            lastJob = runnerSupport.EmptyWorkflowReturn()
        else:
            clobber = lastJob.clobber
        self.lastJob = lastJob
        if not bamFile:
            bamFile = lastJob.data
        self.emailAddress = emailAddress
        self.tempDir = jobList.tempDir
        self.sampleName = sampleName
        splitAndTrim = gatkRunners.SplitNCigarReads(sampleName, bamFile, refGenomeFasta, clobber = clobber, outputDirectory = outputDir)
        self.splitAndTrimCommandIndex = jobList.addJob(splitAndTrim.splitAndTrimCommand)
        splitAndTrimJobID = self.submitCommands(mock)
        self.returnData = runnerSupport.WorkflowReturn(splitAndTrim.bamOut, splitAndTrimJobID, splitAndTrim.clobber)
        
    def submitCommands(self, mock = False):
        import genericRunners
        tempDir = self.tempDir
        splitAndTrim = genericRunners.HoffmanJob(self.lastJob.jobID, self.splitAndTrimCommandIndex, self.sampleName + "splitTrim", tempDir, self.emailAddress, "a", 1, memory = 16, mock = mock)
        return splitAndTrim.jobID
    
class HaplotypeCallerAndVQSR(object):
    
    def __init__(self, jobList, sampleName, refGenomeFasta, jointGenotypingFrozenGVCFDirectory, snpResources, indelResources, bamFile = False, intervals = False, dbSNP = False, emailAddress = False, clobber = None, outputDir = "", lastJob = False, mock = False):
        import houseKeeping
        import gatkRunners
        import runnerSupport
        if not bamFile and not lastJob:
            raise RuntimeError("No bam for processing specified in last job data or arguments.")
        if not lastJob:
            lastJob = runnerSupport.EmptyWorkflowReturn()
        else:
            clobber = lastJob.clobber
        self.lastJob = lastJob
        if not bamFile:
            bamFile = lastJob.data
        self.emailAddress = emailAddress
        self.tempDir = jobList.tempDir
        self.sampleName = sampleName
        haplotypeCaller = gatkRunners.HaplotypeCallerQueue(sampleName, bamFile, refGenomeFasta, intervals, dbSNP = dbSNP, cores = 3, clobber = clobber, outputDirectory = outputDir)
        self.haplotypeCallerCommandIndex = jobList.addJob(haplotypeCaller.haplotypeCallerCommand)
        jointGenotype = gatkRunners.JointGenotype(sampleName, haplotypeCaller.gvcfOut, jointGenotypingFrozenGVCFDirectory, refGenomeFasta, dbSNP, clobber = haplotypeCaller.clobber, outputDirectory = outputDir)
        self.jointGenotypeCommandIndex = jobList.addJob(jointGenotype.jointGenotypeCommand)
        vqsrSNP = gatkRunners.VQSR(sampleName, jointGenotype.vcfOut, snpResources, refGenomeFasta, mode = "SNP", clobber = jointGenotype.clobber, outputDirectory = outputDir)
        self.vqsrSNPAnalysisCommandIndex = jobList.addJob(vqsrSNP.analysisCommand)
        self.vqsrSNPExecuteCommandIndex =jobList.addJob(vqsrSNP.executionCommand)
        vqsrIndel = gatkRunners.VQSR(sampleName, vqsrSNP.vcfOut, indelResources, refGenomeFasta, mode = "INDEL", clobber = vqsrSNP.clobber, outputDirectory = outputDir)
        self.vqsrIndelAnalysisCommandIndex = jobList.addJob(vqsrIndel.analysisCommand)
        self.vqsrIndelExecuteCommandIndex =jobList.addJob(vqsrIndel.executionCommand)
        deleteSNPRecalOnlyVCF = houseKeeping.Delete(vqsrSNP.vcfOut)
        self.deleteSNPRecalOnlyVCFCommandIndex = jobList.addJob(deleteSNPRecalOnlyVCF.deleteCommand)
        vqsrIndelExecuteCommandJobID = self.submitCommands(mock)
        self.returnData = runnerSupport.WorkflowReturn(vqsrIndel.vcfOut, vqsrIndelExecuteCommandJobID, vqsrIndel.clobber)
        
    def submitCommands(self, mock = False):
        import genericRunners
        tempDir = self.tempDir
        haplotypeCaller = genericRunners.HoffmanJob(self.lastJob.jobID, self.haplotypeCallerCommandIndex, self.sampleName + "hapCall", tempDir, self.emailAddress, "a", 1, memory = 20, mock = mock, forceHighP = True)
        jointGenotyper = genericRunners.HoffmanJob(haplotypeCaller.jobID, self.jointGenotypeCommandIndex, self.sampleName + "jointGeno", tempDir, self.emailAddress, "a", 1, memory = 16, mock = mock)
        vqsrSNPAnalysis = genericRunners.HoffmanJob(jointGenotyper.jobID, self.vqsrSNPAnalysisCommandIndex, self.sampleName + "vqsrSNP1", tempDir, self.emailAddress, "a", 1, memory = 16, mock = mock)
        vqsrSNPExecute = genericRunners.HoffmanJob(vqsrSNPAnalysis.jobID, self.vqsrSNPExecuteCommandIndex, self.sampleName + "vqsrSNP2", tempDir, self.emailAddress, "a", 1, memory = 16, mock = mock)
        vqsrIndelAnalysis = genericRunners.HoffmanJob(vqsrSNPExecute.jobID, self.vqsrIndelAnalysisCommandIndex, self.sampleName + "vqsrIndel1", tempDir, self.emailAddress, "a", 1, memory = 16, mock = mock)
        vqsrIndelExecute = genericRunners.HoffmanJob(vqsrIndelAnalysis.jobID, self.vqsrIndelExecuteCommandIndex, self.sampleName + "vqsrIndel2", tempDir, self.emailAddress, "a", 1, memory = 16, mock = mock)
        deleteSNPOnlyRecalVCF = genericRunners.HoffmanJob(vqsrIndelExecute.jobID, self.deleteSNPRecalOnlyVCFCommandIndex, self.sampleName + "delSNPRecalVCF", tempDir, self.emailAddress, "a", 1, mock = mock)
        return vqsrIndelExecute.jobID
    
class HaplotypeCaller(object):
    
    def __init__(self, jobList, sampleName, refGenomeFasta, bamFile = False, intervals = False, dbSNP = False, emailAddress = False, clobber = None, outputDir = "", lastJob = False, mock = False):
        import houseKeeping
        import gatkRunners
        import runnerSupport
        if not bamFile and not lastJob:
            raise RuntimeError("No bam for processing specified in last job data or arguments.")
        if not lastJob:
            lastJob = runnerSupport.EmptyWorkflowReturn()
        else:
            clobber = lastJob.clobber
        self.lastJob = lastJob
        if not bamFile:
            bamFile = lastJob.data
        self.emailAddress = emailAddress
        self.tempDir = jobList.tempDir
        self.sampleName = sampleName
        haplotypeCaller = gatkRunners.HaplotypeCallerQueue(sampleName, bamFile, refGenomeFasta, intervals, dbSNP = dbSNP, cores = 3, clobber = clobber, outputDirectory = outputDir)
        self.haplotypeCallerCommandIndex = jobList.addJob(haplotypeCaller.haplotypeCallerCommand)
        haplotypeCallerJobID = self.submitCommands(mock)
        self.returnData = runnerSupport.WorkflowReturn(haplotypeCaller.gvcfOut, haplotypeCallerJobID, haplotypeCaller.clobber)
        
    def submitCommands(self, mock = False):
        import genericRunners
        tempDir = self.tempDir
        haplotypeCaller = genericRunners.HoffmanJob(self.lastJob.jobID, self.haplotypeCallerCommandIndex, self.sampleName + "hapCall", tempDir, self.emailAddress, "a", 1, memory = 20, mock = mock)
        return haplotypeCaller.jobID
  
class JointGenotype(object):
    
    def __init__(self, jobList, sampleName, refGenomeFasta, jointGenotypingFrozenGVCFDirectory = False, gvcfs = False, intervals = False, dbSNP = False, emailAddress = False, standCallConf = False, clobber = None, outputDir = "", lastJob = False, mock = False):
        import houseKeeping
        import gatkRunners
        import runnerSupport
        if not gvcfs and not lastJob:
            raise RuntimeError("No bam for processing specified in last job data or arguments.")
        if not type(gvcfs) in (bool, list, tuple):
            gvcfs = [gvcfs]
        if not lastJob:
            lastJob = [runnerSupport.EmptyWorkflowReturn()]
        else:
            if not type(lastJob) in (list, tuple):
                lastJob = [lastJob]
            clobberCollector = []
            for job in lastJob:
                clobberCollector.append(job.clobber)
            clobberSet = set(clobberCollector)
            if True in clobberSet:
                clobber = True
            elif False in clobberSet:
                clobber = False
            else:
                clobber = None
        self.dependencies = [item.jobID for item in lastJob]
        if not gvcfs:
            gvcfs = [item.data for item in lastJob]
        self.emailAddress = emailAddress
        self.tempDir = jobList.tempDir
        self.sampleName = sampleName
        jointGenotype = gatkRunners.JointGenotype(sampleName, gvcfs, jointGenotypingFrozenGVCFDirectory, refGenomeFasta, dbSNP, standCallConf = standCallConf, clobber = clobber, outputDirectory = outputDir)
        self.jointGenotypeCommandIndex = jobList.addJob(jointGenotype.jointGenotypeCommand)
        jointGenotyperJobID = self.submitCommands(mock)
        self.returnData = runnerSupport.WorkflowReturn(jointGenotype.vcfOut, jointGenotyperJobID, jointGenotype.clobber)
        
    def submitCommands(self, mock = False):
        import genericRunners
        tempDir = self.tempDir
        jointGenotyper = genericRunners.HoffmanJob(self.dependencies, self.jointGenotypeCommandIndex, self.sampleName + "jointGeno", tempDir, self.emailAddress, "a", 1, memory = 16, mock = mock)
        return jointGenotyper.jobID
    
class JointGenotypeAndVQSR(object):
    
    def __init__(self, jobList, sampleName, refGenomeFasta, jointGenotypingFrozenGVCFDirectory, snpResources, indelResources, gvcfs = False, intervals = False, dbSNP = False, emailAddress = False, clobber = None, outputDir = "", lastJob = False, mock = False):
        import houseKeeping
        import gatkRunners
        import runnerSupport
        if not gvcfs and not lastJob:
            raise RuntimeError("No bam for processing specified in last job data or arguments.")
        if not type(gvcfs) in (bool, list, tuple):
            gvcfs = [gvcfs]
        if not lastJob:
            lastJob = [runnerSupport.EmptyWorkflowReturn()]
        else:
            if not type(lastJob) in (list, tuple):
                lastJob = [lastJob]
            clobberCollector = []
            for job in lastJob:
                clobberCollector.append(job.clobber)
            clobberSet = set(clobberCollector)
            if True in clobberSet:
                clobber = True
            elif False in clobberSet:
                clobber = False
            else:
                clobber = None
        self.dependencies = [item.jobID for item in lastJob]
        if not gvcfs:
            gvcfs = [item.data for item in lastJob]
        self.emailAddress = emailAddress
        self.tempDir = jobList.tempDir
        self.sampleName = sampleName
        jointGenotype = gatkRunners.JointGenotype(sampleName, gvcfs, jointGenotypingFrozenGVCFDirectory, refGenomeFasta, dbSNP, clobber = clobber, outputDirectory = outputDir)
        self.jointGenotypeCommandIndex = jobList.addJob(jointGenotype.jointGenotypeCommand)
        vqsrSNP = gatkRunners.VQSR(sampleName, jointGenotype.vcfOut, snpResources, refGenomeFasta, mode = "SNP", clobber = jointGenotype.clobber, outputDirectory = outputDir)
        self.vqsrSNPAnalysisCommandIndex = jobList.addJob(vqsrSNP.analysisCommand)
        self.vqsrSNPExecuteCommandIndex =jobList.addJob(vqsrSNP.executionCommand)
        vqsrIndel = gatkRunners.VQSR(sampleName, vqsrSNP.vcfOut, indelResources, refGenomeFasta, mode = "INDEL", clobber = vqsrSNP.clobber, outputDirectory = outputDir)
        self.vqsrIndelAnalysisCommandIndex = jobList.addJob(vqsrIndel.analysisCommand)
        self.vqsrIndelExecuteCommandIndex =jobList.addJob(vqsrIndel.executionCommand)
        deleteSNPRecalOnlyVCF = houseKeeping.Delete(vqsrSNP.vcfOut)
        self.deleteSNPRecalOnlyVCFCommandIndex = jobList.addJob(deleteSNPRecalOnlyVCF.deleteCommand)
        vqsrIndelExecuteCommandJobID = self.submitCommands(mock)
        self.returnData = runnerSupport.WorkflowReturn(vqsrIndel.vcfOut, vqsrIndelExecuteCommandJobID, vqsrIndel.clobber)
        
    def submitCommands(self, mock = False):
        import genericRunners
        tempDir = self.tempDir
        jointGenotyper = genericRunners.HoffmanJob(self.dependencies, self.jointGenotypeCommandIndex, self.sampleName + "jointGeno", tempDir, self.emailAddress, "a", 1, memory = 16, mock = mock)
        vqsrSNPAnalysis = genericRunners.HoffmanJob(jointGenotyper.jobID, self.vqsrSNPAnalysisCommandIndex, self.sampleName + "vqsrSNP1", tempDir, self.emailAddress, "a", 1, memory = 16, mock = mock)
        vqsrSNPExecute = genericRunners.HoffmanJob(vqsrSNPAnalysis.jobID, self.vqsrSNPExecuteCommandIndex, self.sampleName + "vqsrSNP2", tempDir, self.emailAddress, "a", 1, memory = 16, mock = mock)
        vqsrIndelAnalysis = genericRunners.HoffmanJob(vqsrSNPExecute.jobID, self.vqsrIndelAnalysisCommandIndex, self.sampleName + "vqsrIndel1", tempDir, self.emailAddress, "a", 1, memory = 16, mock = mock)
        vqsrIndelExecute = genericRunners.HoffmanJob(vqsrIndelAnalysis.jobID, self.vqsrIndelExecuteCommandIndex, self.sampleName + "vqsrIndel2", tempDir, self.emailAddress, "a", 1, memory = 16, mock = mock)
        deleteSNPOnlyRecalVCF = genericRunners.HoffmanJob(vqsrIndelExecute.jobID, self.deleteSNPRecalOnlyVCFCommandIndex, self.sampleName + "delSNPRecalVCF", tempDir, self.emailAddress, "a", 1, mock = mock)
        return vqsrIndelExecute.jobID
    
class Mutect1(object):
    
    def __init__(self, jobList, sampleName, refGenomeFasta, normalBAM = False, tumorBAM = False, intervals = False, dbSNP = False, emailAddress = False, clobber = None, outputDir = "", lastNormalJob = False, lastTumorJob = False, mock = False):
        import houseKeeping
        import gatkRunners
        import programRunners
        import runnerSupport
        if not normalBAM and not lastNormalJob:
            raise RuntimeError("No normal bam for processing specified in last job data or arguments.")
        if not lastNormalJob:
            lastNormalJob = runnerSupport.EmptyWorkflowReturn()
        else:
            clobberNormal = lastNormalJob.clobber
        if not tumorBAM and not lastTumorJob:
            raise RuntimeError("No tumor bam for processing specified in last job data or arguments.")
        if not lastTumorJob:
            lastTumorJob = runnerSupport.EmptyWorkflowReturn()
        else:
            clobberTumor = lastNormalJob.clobber
        if clobber == None:
            if clobberNormal or clobberTumor:
                clobber = True
            elif clobberNormal == False or clobberTumor == False:
                clobber = False
            else:
                clobber = None #not absolutely necessary, since this was the initial condition
        if not normalBAM:
            normalBAM = lastNormalJob.data
        if not tumorBAM:
            tumorBAM = lastTumorJob.data
        self.lastNormalJob = lastNormalJob
        self.lastTumorJob = lastTumorJob
        self.emailAddress = emailAddress
        self.tempDir = jobList.tempDir
        self.sampleName = sampleName
        mutect1 = gatkRunners.Mutect1(sampleName, normalBAM, tumorBAM, refGenomeFasta, intervals, dbSNP, clobber = clobber, outputDirectory = outputDir)
        self.mutect1CommandIndex = jobList.addJob(mutect1.mutect1Command)
        rejectFilter = programRunners.GrepFilter(sampleName, mutect1.vcfOut, mutect1.vcfOut + ".kept", "REJECT", filterOut = True, clobber = mutect1.clobber, outputDirectory = outputDir)
        self.rejectFilterCommandIndex = jobList.addJob(rejectFilter.grepCommand)
        rejectFilterJobID = self.submitCommands(mock)
        self.returnData = runnerSupport.WorkflowReturn(rejectFilter.outputFile, rejectFilterJobID, rejectFilter.clobber)
        
    def submitCommands(self, mock):
        import genericRunners
        tempDir = self.tempDir
        mutect1 = genericRunners.HoffmanJob([self.lastNormalJob.jobID, self.lastTumorJob.jobID], self.mutect1CommandIndex, self.sampleName + "mutect1", tempDir, self.emailAddress, "a", 1, memory = 16, mock = mock)
        rejectFilter = genericRunners.HoffmanJob(mutect1.jobID, self.rejectFilterCommandIndex, self.sampleName + "rejectFilter", tempDir, self.emailAddress, "a", 1, mock = mock)
        return rejectFilter.jobID
    
class Mutect1ToSomaticsData(object):
    
    def __init__(self, jobList, sampleName, refGenomeFasta, normalBAM = False, tumorBAM = False, intervals = False, dbSNP = False, emailAddress = False, clobber = None, outputDir = "", lastNormalJob = False, lastTumorJob = False, minDepth = 10, mock = False):
        import houseKeeping
        import gatkRunners
        import programRunners
        import runnerSupport
        if not normalBAM and not lastNormalJob:
            raise RuntimeError("No normal bam for processing specified in last job data or arguments.")
        if not lastNormalJob:
            lastNormalJob = runnerSupport.EmptyWorkflowReturn()
        else:
            clobberNormal = lastNormalJob.clobber
        if not tumorBAM and not lastTumorJob:
            raise RuntimeError("No tumor bam for processing specified in last job data or arguments.")
        if not lastTumorJob:
            lastTumorJob = runnerSupport.EmptyWorkflowReturn()
        else:
            clobberTumor = lastNormalJob.clobber
        if clobber == None:
            if clobberNormal or clobberTumor:
                clobber = True
            elif clobberNormal == False or clobberTumor == False:
                clobber = False
            else:
                clobber = None #not absolutely necessary, since this was the initial condition
        if not normalBAM:
            normalBAM = lastNormalJob.data
        if not tumorBAM:
            tumorBAM = lastTumorJob.data
        self.lastNormalJob = lastNormalJob
        self.lastTumorJob = lastTumorJob
        self.emailAddress = emailAddress
        self.tempDir = jobList.tempDir
        self.sampleName = sampleName
        mutectAndFilter = Mutect1(jobList, sampleName, refGenomeFasta, normalBAM = normalBAM, tumorBAM = tumorBAM, intervals = intervals, dbSNP = dbSNP, emailAddress = emailAddress, clobber = clobber, outputDir = outputDir, lastNormalJob = lastNormalJob, lastTumorJob = lastTumorJob, mock = mock)
        mutectReader = programRunners.MutectReader(sampleName, mutectAndFilter.returnData.data, minDepth = minDepth, clobber = mutectAndFilter.returnData.clobber, outputDirectory = outputDir)
        self.mutectReaderJobIndex = jobList.addJob(mutectReader.mutectReaderCommand)
        mutectReaderJobID = self.submitCommands(mock, mutectAndFilter)
        self.returnData = runnerSupport.WorkflowReturn(mutectReader.acceptedOut, mutectReaderJobID, mutectReader.clobber)
        
    def submitCommands(self, mock, mutectAndFilter):
        import genericRunners
        tempDir = self.tempDir
        mutectReader = genericRunners.HoffmanJob([mutectAndFilter.returnData.jobID], self.mutectReaderJobIndex, self.sampleName + "MuRead", tempDir, self.emailAddress, "a", 1, mock = mock)
        return mutectReader.jobID
    
class HaplotypeCallerSomaticsData(object):
    
    def __init__(self, jobList, sampleName, refGenomeFasta, jointGenotypingFrozenGVCFDirectory, snpResources, indelResources, normalBAM = False, tumorBAM = False, intervals = False, dbSNP = False, emailAddress = False, clobber = None, outputDir = "", lastNormalJob = False, lastTumorJob = False, maxPValue = 0.05, minDepth = 10, mock = False):
        import houseKeeping
        import gatkRunners
        import programRunners
        import runnerSupport
        if not normalBAM and not lastNormalJob:
            raise RuntimeError("No normal bam for processing specified in last job data or arguments.")
        if not lastNormalJob:
            lastNormalJob = runnerSupport.EmptyWorkflowReturn()
        else:
            clobberNormal = lastNormalJob.clobber
        if not tumorBAM and not lastTumorJob:
            raise RuntimeError("No tumor bam for processing specified in last job data or arguments.")
        if not lastTumorJob:
            lastTumorJob = runnerSupport.EmptyWorkflowReturn()
        else:
            clobberTumor = lastNormalJob.clobber
        if clobber == None:
            if clobberNormal or clobberTumor:
                clobber = True
            elif clobberNormal == False or clobberTumor == False:
                clobber = False
            else:
                clobber = None #not absolutely necessary, since this was the initial condition
        if not normalBAM:
            normalBAM = lastNormalJob.data
        if not tumorBAM:
            tumorBAM = lastTumorJob.data
        tumorBaseName = runnerSupport.baseFileName(tumorBAM)
        normalBaseName = runnerSupport.baseFileName(normalBAM)
        self.lastNormalJob = lastNormalJob
        self.lastTumorJob = lastTumorJob
        self.emailAddress = emailAddress
        self.tempDir = jobList.tempDir
        self.sampleName = sampleName
        normalHaplotypeCaller = HaplotypeCaller(jobList, normalBaseName, refGenomeFasta, intervals = intervals, dbSNP = dbSNP, emailAddress = emailAddress, clobber = clobber, outputDir = outputDir, lastJob = lastNormalJob, mock = mock)
        tumorHaplotypeCaller = HaplotypeCaller(jobList, tumorBaseName, refGenomeFasta, intervals = intervals, dbSNP = dbSNP, emailAddress = emailAddress, clobber = normalHaplotypeCaller.returnData.clobber, outputDir = outputDir, lastJob = lastTumorJob, mock = mock)
        jointGenotyping = JointGenotypeAndVQSR(jobList, sampleName, refGenomeFasta, jointGenotypingFrozenGVCFDirectory, snpResources, indelResources, gvcfs = [normalHaplotypeCaller.returnData.data, tumorHaplotypeCaller.returnData.data], intervals = intervals, dbSNP = dbSNP, emailAddress = emailAddress, clobber = tumorHaplotypeCaller.returnData.clobber, outputDir = outputDir, lastJob = [normalHaplotypeCaller.returnData, tumorHaplotypeCaller.returnData], mock = mock)
        vcfReader = programRunners.VCFReader(sampleName, jointGenotyping.returnData.data, tumorBaseName, normalBaseName, maxPValue, minDepth, jointGenotyping.returnData.clobber, outputDirectory = outputDir)
        self.vcfReaderJobIndex = jobList.addJob(vcfReader.vcfReaderCommand)
        vcfReaderJobID = self.submitCommands(mock, jointGenotyping)
        self.returnData = runnerSupport.WorkflowReturn(vcfReader.acceptedOut, vcfReaderJobID, vcfReader.clobber)
        
class VCFReader(object):
    
    def __init__(self, jobList, normalSampleName, tumorSampleName, vcf = False, emailAddress = False, clobber = None, outputDir = "", lastJob = False, maxPValue = 0.05, minDepth = 10, mock = False):
        import houseKeeping
        import gatkRunners
        import programRunners
        import runnerSupport
        if clobber == None:
            clobber = lastJob.clobber
        else:
            clobber = clobber
        if lastJob:
            self.vcf = lastJob.data
        else:
            self.vcf = vcf
            lastJob = runnerSupport.EmptyWorkflowReturn()
        if not self.vcf:
            raise RuntimeError("Error: No VCF to analyze was given.")
        self.lastJob = lastJob
        self.emailAddress = emailAddress
        self.tempDir = jobList.tempDir
        self.sampleName = tumorSampleName
        vcfReader = programRunners.VCFReader(tumorSampleName, self.vcf, tumorSampleName, normalSampleName, maxPValue, minDepth, clobber, outputDirectory = outputDir)
        self.vcfReaderJobIndex = jobList.addJob(vcfReader.vcfReaderCommand)
        vcfReaderJobID = self.submitCommands(mock)
        self.returnData = runnerSupport.WorkflowReturn(vcfReader.acceptedOut, vcfReaderJobID, vcfReader.clobber)
        
    def submitCommands(self, mock):
        import genericRunners
        tempDir = self.tempDir
        vcfReader = genericRunners.HoffmanJob(self.lastJob.jobID, self.vcfReaderJobIndex, self.sampleName + "VCFRead", tempDir, self.emailAddress, "a", 1, mock = mock)
        return vcfReader.jobID
    
class HaplotypeCallerRNAAnalysis(object):
    
    def __init__(self, jobList, sampleName, refGenomeFasta, rnaBAM = False, tumorBAM = False, intervals = False, dbSNP = False, emailAddress = False, clobber = None, outputDir = "", lastRNAJob = False, lastTumorJob = False, mock = False):
        import houseKeeping
        import gatkRunners
        import programRunners
        import runnerSupport
        sampleName = sampleName + "EXP" #tagging these as expression to avoid collisions with non-expression jobs
        if not rnaBAM and not lastRNAJob:
            raise RuntimeError("No RNA bam for processing specified in last job data or arguments.")
        if not lastRNAJob:
            lastRNAJob = runnerSupport.EmptyWorkflowReturn()
        else:
            clobberRNA = lastRNAJob.clobber
        if not tumorBAM and not lastTumorJob:
            raise RuntimeError("No tumor bam for processing specified in last job data or arguments.")
        if not lastTumorJob:
            lastTumorJob = runnerSupport.EmptyWorkflowReturn()
        else:
            clobberTumor = lastRNAJob.clobber
        if clobber == None:
            if clobberRNA or clobberTumor:
                clobber = True
            elif clobberRNA == False or clobberTumor == False:
                clobber = False
            else:
                clobber = None #not absolutely necessary, since this was the initial condition
        if not rnaBAM:
            rnaBAM = lastRNAJob.data
        if not tumorBAM:
            tumorBAM = lastTumorJob.data
        tumorBaseName = runnerSupport.baseFileName(tumorBAM)
        rnaBaseName = runnerSupport.baseFileName(rnaBAM)
        self.lastRNAJob = lastRNAJob
        self.lastTumorJob = lastTumorJob
        self.emailAddress = emailAddress
        self.tempDir = jobList.tempDir
        self.sampleName = sampleName
        tumorRNAExomeHaplotypeCall = gatkRunners.HaplotypeCallerQueue(self.sampleName, [tumorBAM, rnaBAM], refGenomeFasta, bedFile = intervals, emitRefConfidence = False, standCallConf = 20, variantIndexType = False, variantIndexParameter = False, dbSNP = dbSNP, intervalPadding = 100, allowPotentiallyMisencodedQualityScores = True, cores = 3, clobber = clobber, outputDirectory = outputDir)
        self.haplotypeCallerCommandIndex = jobList.addJob(tumorRNAExomeHaplotypeCall.haplotypeCallerCommand)
        haplotypeCallerJobID = self.submitCommands(mock)
        self.returnData = runnerSupport.WorkflowReturn(tumorRNAExomeHaplotypeCall.vcfOut, haplotypeCallerJobID, tumorRNAExomeHaplotypeCall.clobber)
        
    def submitCommands(self, mock):
        import genericRunners
        tempDir = self.tempDir
        haplotypeCaller = genericRunners.HoffmanJob([self.lastRNAJob.jobID, self.lastTumorJob.jobID], self.haplotypeCallerCommandIndex, self.sampleName + "callTumorRNA", tempDir, self.emailAddress, "a", 3, memory = 20, mock = mock)
        return haplotypeCaller.jobID
    
class VarScan(object):
    
    def __init__(self, jobList, sampleName, refGenomeFasta, normalBAM = False, tumorBAM = False, emailAddress = False, clobber = None, outputDir = "", lastNormalJob = False, lastTumorJob = False, mock = False):
        import houseKeeping
        import gatkRunners
        import programRunners
        import runnerSupport
        if not normalBAM and not lastNormalJob:
            raise RuntimeError("No normal bam for processing specified in last job data or arguments.")
        if not lastNormalJob:
            lastNormalJob = runnerSupport.EmptyWorkflowReturn()
        else:
            clobberNormal = lastNormalJob.clobber
        if not tumorBAM and not lastTumorJob:
            raise RuntimeError("No tumor bam for processing specified in last job data or arguments.")
        if not lastTumorJob:
            lastTumorJob = runnerSupport.EmptyWorkflowReturn()
        else:
            clobberTumor = lastNormalJob.clobber
        if clobber == None:
            if clobberNormal or clobberTumor:
                clobber = True
            elif clobberNormal == False or clobberTumor == False:
                clobber = False
            else:
                clobber = None #not absolutely necessary, since this was the initial condition
        if not normalBAM:
            normalBAM = lastNormalJob.data
        if not tumorBAM:
            tumorBAM = lastTumorJob.data
        self.lastNormalJob = lastNormalJob
        self.lastTumorJob = lastTumorJob
        self.emailAddress = emailAddress
        self.tempDir = jobList.tempDir
        self.sampleName = sampleName
        mPileupNormal = programRunners.MPileup(sampleName + "norm", normalBAM, refGenomeFasta, bgzip = False, clobber = clobber, outputDirectory = outputDir)
        self.mPileupNormalCommandIndex = jobList.addJob(mPileupNormal.mPileupCommand)
        mPileupTumor = programRunners.MPileup(sampleName + "tumor", tumorBAM, refGenomeFasta, bgzip = False, clobber = mPileupNormal.clobber, outputDirectory = outputDir)
        self.mPileupTumorCommandIndex = jobList.addJob(mPileupTumor.mPileupCommand)
        varScan = programRunners.Varscan(sampleName, mPileupNormal.mPileupOut, mPileupTumor.mPileupOut, clobber = mPileupTumor.clobber, outputDirectory = outputDir)
        self.varScanSomaticCommandIndex = jobList.addJob(varScan.varscanSomaticCommand)
        self.varScanProcessSNPCommandIndex = jobList.addJob(varScan.varscanProcessSNPCommand)
        self.varScanProcessIndelCommandIndex = jobList.addJob(varScan.varscanProcessIndelCommand)
        snpJobID, indelJobID = self.submitCommands(mock)
        self.returnData = {"snp":runnerSupport.WorkflowReturn(varScan.hcSomaticSNPOut, snpJobID, varScan.clobber),
                      "indel":runnerSupport.WorkflowReturn(varScan.hcSomaticIndelOut, indelJobID, varScan.clobber)}
        
    def submitCommands(self, mock):
        import genericRunners
        tempDir = self.tempDir
        mPileupNormal = genericRunners.HoffmanJob([self.lastNormalJob.jobID], self.mPileupNormalCommandIndex, self.sampleName + "normPile", tempDir, self.emailAddress, "a", 1, mock = mock)
        mPileupTumor = genericRunners.HoffmanJob([self.lastTumorJob.jobID], self.mPileupTumorCommandIndex, self.sampleName + "tumorPile", tempDir, self.emailAddress, "a", 1, mock = mock)
        varScanSomatic = genericRunners.HoffmanJob([mPileupNormal.jobID, mPileupTumor.jobID], self.varScanSomaticCommandIndex, self.sampleName + "varScanSom", tempDir, self.emailAddress, "a", 1, mock = mock)
        varScanSNP = genericRunners.HoffmanJob([varScanSomatic.jobID], self.varScanProcessSNPCommandIndex, self.sampleName + "varScanSNP", tempDir, self.emailAddress, "a", 1, mock = mock)
        varScanIndel = genericRunners.HoffmanJob([varScanSomatic.jobID], self.varScanProcessIndelCommandIndex, self.sampleName  +"varScanIndel", tempDir, self.emailAddress, "a", 1, mock = mock)
        return [varScanSNP.jobID, varScanIndel.jobID]
    
class VarScanSomaticsData(object):
    
    def __init__(self, jobList, sampleName, refGenomeFasta, normalBAM = False, tumorBAM = False, emailAddress = False, clobber = None, outputDir = "", lastNormalJob = False, lastTumorJob = False, maxWarnings = 10, minDepth = 10, mock = False):
        import houseKeeping
        import gatkRunners
        import programRunners
        import runnerSupport
        if not normalBAM and not lastNormalJob:
            raise RuntimeError("No normal bam for processing specified in last job data or arguments.")
        if not lastNormalJob:
            lastNormalJob = runnerSupport.EmptyWorkflowReturn()
        else:
            clobberNormal = lastNormalJob.clobber
        if not tumorBAM and not lastTumorJob:
            raise RuntimeError("No tumor bam for processing specified in last job data or arguments.")
        if not lastTumorJob:
            lastTumorJob = runnerSupport.EmptyWorkflowReturn()
        else:
            clobberTumor = lastNormalJob.clobber
        if clobber == None:
            if clobberNormal or clobberTumor:
                clobber = True
            elif clobberNormal == False or clobberTumor == False:
                clobber = False
            else:
                clobber = None #not absolutely necessary, since this was the initial condition
        if not normalBAM:
            normalBAM = lastNormalJob.data
        if not tumorBAM:
            tumorBAM = lastTumorJob.data
        self.lastNormalJob = lastNormalJob
        self.lastTumorJob = lastTumorJob
        self.emailAddress = emailAddress
        self.tempDir = jobList.tempDir
        self.sampleName = sampleName
        varScan = VarScan(jobList, sampleName, refGenomeFasta, normalBAM = normalBAM, tumorBAM = tumorBAM, emailAddress = emailAddress, clobber = clobber, outputDir = outputDir, lastNormalJob = lastNormalJob, lastTumorJob = lastTumorJob, mock = mock)
        snpPositionPuller = programRunners.VarScanPositions(sampleName, varScan.returnData["snp"].data, clobber = varScan.returnData["snp"].clobber, outputDirectory = outputDir)
        self.snpPositionPullerCommandIndex = jobList.addJob(snpPositionPuller.positionPullerCommand)
        indelPositionPuller = programRunners.VarScanPositions(sampleName, varScan.returnData["indel"].data, clobber = snpPositionPuller.clobber, outputDirectory = outputDir)
        self.indelPositionPullerCommandIndex = jobList.addJob(indelPositionPuller.positionPullerCommand)
        snpBAMReadCount = programRunners.BAMReadCount(sampleName, tumorBAM, snpPositionPuller.positionsOut, refGenomeFasta, maxWarnings = maxWarnings, clobber = indelPositionPuller.clobber, outputDirectory = outputDir)
        self.snpBAMReadCountCommandIndex = jobList.addJob(snpBAMReadCount.bamReadcountCommand)
        indelBAMReadCount = programRunners.BAMReadCount(sampleName, tumorBAM, indelPositionPuller.positionsOut, refGenomeFasta, maxWarnings = maxWarnings, clobber = snpBAMReadCount.clobber, outputDirectory = outputDir)
        self.indelBAMReadCountCommandIndex = jobList.addJob(indelBAMReadCount.bamReadcountCommand)
        snpFPFilter = programRunners.VarScanFPFilter(sampleName, varScan.returnData["snp"].data, snpBAMReadCount.readcountOut, clobber = indelBAMReadCount.clobber, outputDirectory = outputDir)
        self.snpFPFilterCommandIndex = jobList.addJob(snpFPFilter.fpFilterCommand)
        indelFPFilter = programRunners.VarScanFPFilter(sampleName, varScan.returnData["indel"].data, indelBAMReadCount.readcountOut, clobber = snpFPFilter.clobber, outputDirectory = outputDir)
        self.indelFPFilterCommandIndex = jobList.addJob(indelFPFilter.fpFilterCommand)
        snpVarScanReader = programRunners.VarScanReader(sampleName, snpFPFilter.passOut, minDepth = minDepth, clobber = clobber, outputDirectory = outputDir)
        self.snpVarScanReaderCommandIndex = jobList.addJob(snpVarScanReader.varScanReaderCommand)
        indelVarScanReader = programRunners.VarScanReader(sampleName, indelFPFilter.passOut, minDepth = minDepth, clobber = clobber, outputDirectory = outputDir)
        self.indelVarScanReaderCommandIndex = jobList.addJob(indelVarScanReader.varScanReaderCommand)
        snpJobID, indelJobID = self.submitCommands(mock, varScan)
        self.returnData = {"snp":runnerSupport.WorkflowReturn(snpVarScanReader.acceptedOut, snpJobID, snpVarScanReader.clobber),
                           "indel":runnerSupport.WorkflowReturn(indelVarScanReader.acceptedOut, indelJobID, indelVarScanReader.clobber)}
        
    def submitCommands(self, mock, varScan):
        import genericRunners
        tempDir = self.tempDir
        snpPositionPuller = genericRunners.HoffmanJob([varScan.returnData["snp"].jobID], self.snpPositionPullerCommandIndex, self.sampleName + "snpPos", tempDir, self.emailAddress, "a", 1, mock = mock)
        indelPositionPuller = genericRunners.HoffmanJob([varScan.returnData["indel"].jobID], self.indelPositionPullerCommandIndex, self.sampleName + "indelPos", tempDir, self.emailAddress, "a", 1, mock = mock)
        snpBAMReadCount = genericRunners.HoffmanJob([snpPositionPuller.jobID], self.snpBAMReadCountCommandIndex, self.sampleName + "snpBAMRead", tempDir, self.emailAddress, "a", 1, mock = mock)
        indelBAMReadCount = genericRunners.HoffmanJob([indelPositionPuller.jobID], self.indelBAMReadCountCommandIndex, self.sampleName + "indelBAMRead", tempDir, self.emailAddress, "a", 1, mock = mock)
        snpFPFilter = genericRunners.HoffmanJob([snpBAMReadCount.jobID], self.snpFPFilterCommandIndex, self.sampleName + "snpFPF", tempDir, self.emailAddress, "a", 1, mock = mock)
        indelFPFilter = genericRunners.HoffmanJob([indelBAMReadCount.jobID], self.indelFPFilterCommandIndex, self.sampleName + "indelFPF", tempDir, self.emailAddress, "a", 1, mock = mock)
        snpVarScanReader = genericRunners.HoffmanJob([snpFPFilter.jobID], self.snpVarScanReaderCommandIndex, self.sampleName + "snpVSR", tempDir, self.emailAddress, "a", 1, mock = mock)
        indelVarScanReader = genericRunners.HoffmanJob([indelFPFilter.jobID], self.indelVarScanReaderCommandIndex, self.sampleName + "indelVSR", tempDir, self.emailAddress, "a", 1, mock = mock)
        return [snpVarScanReader.jobID, indelVarScanReader.jobID]
    
class CombineVarScanHaplotypeCallerMutectSomatics(object):
    
    def __init__(self, jobList, sampleName, varScanSomaticIndels = False, varScanSomaticSNPs = False, haplotypeCallerSomatics = False, mutect1Somatics = False, emailAddress = False, clobber = None, outputDir = "", varScanSNPSomaticsJob = False, varScanIndelSomaticsJob = False, varScanCombinedSomaticsJob = False, haplotypeCallerSomaticsJob = False, mutect1SomaticsJob = False, mock = False):
        import houseKeeping
        import gatkRunners
        import programRunners
        import runnerSupport
        if not ((varScanSomaticIndels and varScanSomaticSNPs) or (varScanIndelSomaticsJob and varScanSNPSomaticsJob) or varScanCombinedSomaticsJob):
            raise RuntimeError("There appears to be missing data for VarScan")
        if not (haplotypeCallerSomatics or haplotypeCallerSomaticsJob):
            raise RuntimeError("Missing haplotype caller data.")
        if not (mutect1Somatics or mutect1SomaticsJob):
            raise RuntimeError("Missing mutect1 data")
        if varScanCombinedSomaticsJob:
            varScanIndelSomaticsJob = varScanCombinedSomaticsJob["indel"]
            varScanSNPSomaticsJob = varScanCombinedSomaticsJob["snp"]
        if not varScanIndelSomaticsJob:
            varScanIndelSomaticsJob = runnerSupport.EmptyWorkflowReturn()
        if not varScanSNPSomaticsJob:
            varScanSNPSomaticsJob = runnerSupport.EmptyWorkflowReturn()
        if not haplotypeCallerSomaticsJob:
            haplotypeCallerSomaticsJob = runnerSupport.EmptyWorkflowReturn()
        if not mutect1SomaticsJob:
            mutect1SomaticsJob = runnerSupport.EmptyWorkflowReturn()
        clobberSet = set([clobber] + [item.clobber for item in [haplotypeCallerSomaticsJob, mutect1SomaticsJob, varScanIndelSomaticsJob, varScanSNPSomaticsJob]])
        if True in clobberSet:
            clobber = True
        elif False  in clobberSet:
            clobber = False
        elif None in clobberSet:
            clobber = None
        if not varScanSomaticIndels:
            varScanSomaticIndels = "varscan,indel,%s" %varScanIndelSomaticsJob.data
        if not varScanSomaticSNPs:
            varScanSomaticSNPs = "varscan,snp,%s" %varScanSNPSomaticsJob.data
        if not haplotypeCallerSomatics:
            haplotypeCallerSomatics = "haplotypecaller,%s" %haplotypeCallerSomaticsJob.data
        if not mutect1Somatics:
            mutect1Somatics = "mutect1,%s" %mutect1SomaticsJob.data
        self.varScanIndelSomaticsJob = varScanIndelSomaticsJob
        self.varScanSNPSomaticsJob = varScanSNPSomaticsJob
        self.haplotypeCallerSomaticsJob = haplotypeCallerSomaticsJob
        self.mutect1SomaticsJob = mutect1SomaticsJob
        self.emailAddress = emailAddress
        self.tempDir = jobList.tempDir
        self.sampleName = sampleName
        variantCombine = programRunners.VariantCombine(sampleName, [varScanSomaticIndels, varScanSomaticSNPs, haplotypeCallerSomatics, mutect1Somatics], minHits = 2, maxHits = False, clobber = clobber, outputDirectory = outputDir)
        self.variantCombineCommandIndex = jobList.addJob(variantCombine.variantCombineCommand)
        variantCombineJobID = self.submitCommands(mock)
        self.returnData = runnerSupport.WorkflowReturn(variantCombine.acceptedSomaticVariants, variantCombineJobID, variantCombine.clobber)
        
    def submitCommands(self, mock):
        import genericRunners
        tempDir = self.tempDir
        variantCombine = genericRunners.HoffmanJob([self.varScanIndelSomaticsJob.jobID, self.varScanSNPSomaticsJob.jobID, self.haplotypeCallerSomaticsJob.jobID, self.mutect1SomaticsJob.jobID], self.variantCombineCommandIndex, self.sampleName + "Combine", tempDir, self.emailAddress, "a", 1, mock = mock)
        return variantCombine.jobID
    
class MPileup(object):
    def __init__(self, jobList, sampleName, bamFile = False, refGenomeFasta = False, gzip = False, bgzip = False, emailAddress = False, clobber = None, outputDir = "", makeBAMJob = False, mock = False):
        import houseKeeping
        import gatkRunners
        import programRunners
        import runnerSupport
        if not (makeBAMJob or bamFile):
            raise RuntimeError("No BAM file passed through job info or the explicit argument.  Nothing to work on.")
        if bamFile:
            self.bamFile = bamFile
        else:
            self.bamFile = makeBAMJob.data
        if not makeBAMJob:
            self.makeBAMJob = runnerSupport.EmptyWorkflowReturn()
        else:
            self.makeBAMJob = makeBAMJob
        if makeBAMJob.clobber:
            clobber = makeBAMJob.clobber
        self.emailAddress = emailAddress
        self.tempDir = jobList.tempDir
        self.sampleName = sampleName
        mPileup = programRunners.MPileup(sampleName, self.bamFile, refGenomeFasta, gzip = gzip, bgzip = bgzip, clobber = clobber, outputDirectory = outputDir)
        self.mPileupCommandIndex = jobList.addJob(mPileup.mPileupCommand)
        mPileupJobID = self.submitCommands(mock)
        self.returnData = runnerSupport.WorkflowReturn(mPileup.mPileupOut, mPileupJobID, mPileup.clobber)
        
    def submitCommands(self, mock):
        import genericRunners
        tempDir = self.tempDir
        mPileup = genericRunners.HoffmanJob([self.makeBAMJob.jobID], self.mPileupCommandIndex, self.sampleName + "mPileup", tempDir, self.emailAddress, "a", 1, mock = mock)
        return mPileup.jobID
    
    
class GetRNASupportVCF(object):
    
    def __init__(self, jobList, sampleName, rnaSampleName, somaticVariantPickle = False, rnaVariantVCF = False, minDifference = 10, emailAddress = False, clobber = None, outputDir = "", combineSomaticsJob = False, rnaVariantCallJob = False, mock = False):
        import houseKeeping
        import gatkRunners
        import programRunners
        import runnerSupport
        self.minDifference = minDifference
        if not (somaticVariantPickle or combineSomaticsJob):
            raise RuntimeError("No combined somatics job data or combined somatics file passed.")
        if not (rnaVariantCallJob or rnaVariantVCF):
            raise RuntimeError("No RNA variant call job data or RNA variant call VCF passed.")
        if rnaVariantVCF:
            self.rnaVariantVCF = rnaVariantVCF
        else:
            self.rnaVariantVCF = rnaVariantCallJob.data
        if somaticVariantPickle:
            self.somaticVariantPickle = somaticVariantPickle
        else:
            self.somaticVariantPickle = combineSomaticsJob.data
        if not combineSomaticsJob:
            combineSomaticsJob = runnerSupport.EmptyWorkflowReturn()
        if not rnaVariantCallJob:
            rnaVariantCallJob = runnerSupport.EmptyWorkflowReturn()
        self.combineSomaticsJob = combineSomaticsJob
        self.rnaVariantCallJob = rnaVariantCallJob
        clobberSet = set([clobber] + [item.clobber for item in [combineSomaticsJob, rnaVariantCallJob]])
        if True in clobberSet:
            clobber = True
        elif False  in clobberSet:
            clobber = False
        elif None in clobberSet:
            clobber = None
        self.emailAddress = emailAddress
        self.tempDir = jobList.tempDir
        self.sampleName = sampleName
        getRNASupportVCF = programRunners.GetRNASupportVCF(sampleName, rnaSampleName, self.somaticVariantPickle, self.rnaVariantVCF, minDifference, clobber = clobber, outputDirectory = outputDir)
        self.getRNASupportCommandIndex = jobList.addJob(getRNASupportVCF.getRNASupportCommand)
        getRNASupportJobID = self.submitCommands(mock)
        self.returnData = runnerSupport.WorkflowReturn(getRNASupportVCF.rnaSupportOut, getRNASupportJobID, getRNASupportVCF.clobber)
        
    def submitCommands(self, mock):
        import genericRunners
        tempDir = self.tempDir
        getRNASupport = genericRunners.HoffmanJob([self.rnaVariantCallJob.jobID, self.combineSomaticsJob.jobID], self.getRNASupportCommandIndex, self.sampleName + "AddRNA", tempDir, self.emailAddress, "a", 1, mock = mock)
        return getRNASupport.jobID
    
class GetRNASupportMPileup(object):
    
    def __init__(self, jobList, sampleName, somaticVariantPickle = False, rnaMPileup = False, minDifference = 10, emailAddress = False, clobber = None, outputDir = "", combineSomaticsJob = False, rnaMPileupJob = False, mock = False):
        import houseKeeping
        import gatkRunners
        import programRunners
        import runnerSupport
        self.minDifference = minDifference
        if not (somaticVariantPickle or combineSomaticsJob):
            raise RuntimeError("No combined somatics job data or combined somatics file passed.")
        if not (rnaMPileupJob or rnaMPileup):
            raise RuntimeError("No RNA variant call job data or RNA variant call VCF passed.")
        if rnaMPileup:
            self.rnaMPileup = rnaMPileup
        else:
            self.rnaMPileup = rnaMPileupJob.data
        if somaticVariantPickle:
            self.somaticVariantPickle = somaticVariantPickle
        else:
            self.somaticVariantPickle = combineSomaticsJob.data
        if not combineSomaticsJob:
            combineSomaticsJob = runnerSupport.EmptyWorkflowReturn()
        if not rnaMPileupJob:
            rnaMPileupJob = runnerSupport.EmptyWorkflowReturn()
        self.combineSomaticsJob = combineSomaticsJob
        self.rnaMPileupJob = rnaMPileupJob
        clobberSet = set([clobber] + [item.clobber for item in [combineSomaticsJob, rnaMPileupJob]])
        if True in clobberSet:
            clobber = True
        elif False  in clobberSet:
            clobber = False
        elif None in clobberSet:
            clobber = None
        self.emailAddress = emailAddress
        self.tempDir = jobList.tempDir
        self.sampleName = sampleName
        getRNASupportMPileup = programRunners.GetRNASupportMPileup(sampleName, self.somaticVariantPickle, self.rnaMPileup, minDifference = minDifference, outputFormat = "pickle", clobber = clobber, outputDirectory = outputDir)
        self.getRNASupportCommandIndex = jobList.addJob(getRNASupportMPileup.getRNASupportCommand)
        getRNASupportJobID = self.submitCommands(mock)
        self.returnData = runnerSupport.WorkflowReturn(getRNASupportMPileup.rnaSupportOut, getRNASupportJobID, getRNASupportMPileup.clobber)
        
    def submitCommands(self, mock):
        import genericRunners
        tempDir = self.tempDir
        getRNASupport = genericRunners.HoffmanJob([self.rnaMPileupJob.jobID, self.combineSomaticsJob.jobID], self.getRNASupportCommandIndex, self.sampleName + "AddRNA", tempDir, self.emailAddress, "a", 1, mock = mock)
        return getRNASupport.jobID
    
class CombineTandemVariants(object):
    
    def __init__(self, jobList, sampleName, combinedVariantPickle = False, maxFusionLength = 0, maxDifferencePercent = 10, emailAddress = False, clobber = None, outputDir = "", lastVariantJob = False, mock = False):
        import houseKeeping
        import gatkRunners
        import programRunners
        import runnerSupport
        if not (combinedVariantPickle or lastVariantJob):
            raise RuntimeError("No combined somatics job data or combined somatics file passed.")
        if combinedVariantPickle:
            self.combinedVariantPickle = combinedVariantPickle
        else:
            self.combinedVariantPickle = lastVariantJob.data
        if not lastVariantJob:
            lastVariantJob = runnerSupport.EmptyWorkflowReturn()
        self.lastVariantJob = lastVariantJob
        clobberSet = set([clobber] + [item.clobber for item in [lastVariantJob]])
        if True in clobberSet:
            clobber = True
        elif False  in clobberSet:
            clobber = False
        elif None in clobberSet:
            clobber = None
        self.emailAddress = emailAddress
        self.tempDir = jobList.tempDir
        self.sampleName = sampleName
        fuseVariants = programRunners.TandemVariantCombine(sampleName, self.combinedVariantPickle, outputFormat = "pickle", maxFusionLength = maxFusionLength, maxDifferencePercent = maxDifferencePercent, clobber = clobber, outputDirectory = outputDir)
        self.fusionCommandIndex = jobList.addJob(fuseVariants.fusionCommand)
        fuseVariantsJobID = self.submitCommands(mock)
        self.returnData = runnerSupport.WorkflowReturn(fuseVariants.fusionOut, fuseVariantsJobID, fuseVariants.clobber)
        
    def submitCommands(self, mock):
        import genericRunners
        tempDir = self.tempDir
        variantFusion = genericRunners.HoffmanJob([self.lastVariantJob.jobID], self.fusionCommandIndex, self.sampleName + "fuseVariants", tempDir, self.emailAddress, "a", 1, mock = mock)
        return variantFusion.jobID
    
class GetPeptides(object):
    
    def __init__(self, jobList, sampleName, proteinFasta, fusedVariantPickle = False, sizeList = [8, 9, 10], emailAddress = False, clobber = None, outputDir = "", oncotatorJob = False, mock = False):
        import houseKeeping
        import gatkRunners
        import programRunners
        import runnerSupport
        if not (fusedVariantPickle or oncotatorJob):
            raise RuntimeError("No fused somatics job data or fused somatics file passed.")
        if fusedVariantPickle:
            self.fusedVariantPickle = fusedVariantPickle
        else:
            self.fusedVariantPickle = oncotatorJob.data
        if not oncotatorJob:
            oncotatorJob = runnerSupport.EmptyWorkflowReturn()
        self.oncotatorJob = oncotatorJob
        clobberSet = set([clobber] + [item.clobber for item in [oncotatorJob]])
        if True in clobberSet:
            clobber = True
        elif False  in clobberSet:
            clobber = False
        elif None in clobberSet:
            clobber = None
        self.emailAddress = emailAddress
        self.tempDir = jobList.tempDir
        self.sampleName = sampleName
        getPeptides = programRunners.GetPeptides(sampleName, self.fusedVariantPickle, proteinFasta, sizeList, clobber = clobber, outputDirectory = outputDir)
        self.getPeptidesCommandIndex = jobList.addJob(getPeptides.getPeptidesCommand)
        getPeptidesJobID = self.submitCommands(mock)
        self.returnData = runnerSupport.WorkflowReturn(getPeptides.peptideDataOut, getPeptidesJobID, getPeptides.clobber)
        
    def submitCommands(self, mock):
        import genericRunners
        tempDir = self.tempDir
        getPeptides = genericRunners.HoffmanJob([self.oncotatorJob.jobID], self.getPeptidesCommandIndex, self.sampleName + "peptides", tempDir, self.emailAddress, "a", 1, mock = mock)
        return getPeptides.jobID
    
class SickleTrim(object):
    
    def __init__(self, jobList, sampleName, forwardFastqIn, reverseFastqIn = False, pairedEndMode = True, qualityScoreFormat = "sanger", minWindowQualityAverage = 20, minLengthAfterTrimming = 20, gzippedOutput = True, emailAddress = False, clobber = False, outputDirectory = "", lastForwardJob = False, lastReverseJob = False, mock = False):
        import houseKeeping
        import gatkRunners
        import programRunners
        import runnerSupport
        self.sampleName = sampleName
        self.gzippedOut = gzippedOutput
        self.pairedEndMode = pairedEndMode
        self.qualityScoreFormat = qualityScoreFormat
        self.minWindowQualityAverage = minWindowQualityAverage
        self.minLengthAfterTrimming = minLengthAfterTrimming
        if not (lastForwardJob or forwardFastqIn):
            raise RuntimeError("No forward fastq specified explicitly or from a previous job")
        if pairedEndMode:
            if not (lastReverseJob or reverseFastqIn):
                raise RuntimeError("No reverse fastq specified explicitly or from a previous job and job is running in paired-end mode.")
        if forwardFastqIn:
            self.forwardFastqIn = forwardFastqIn
        else:
            self.forwardFastqIn = lastForwardJob.data
        if pairedEndMode:
            if reverseFastqIn:
                self.reverseFastqIn = reverseFastqIn
            else:
                self.reverseFastqIn = lastReverseJob.data
        if not lastForwardJob:
            self.lastForwardJob = runnerSupport.EmptyWorkflowReturn()
        if not lastReverseJob:
            self.lastReverseJob = runnerSupport.EmptyWorkflowReturn()
        clobberSet = set([clobber] + [item.clobber for item in [self.lastForwardJob, self.lastReverseJob]])
        if True in clobberSet:
            clobber = True
        elif False  in clobberSet:
            clobber = False
        elif None in clobberSet:
            clobber = None
        self.emailAddress = emailAddress
        self.tempDir = jobList.tempDir
        self.sampleName = sampleName
        sickleTrim = programRunners.Sickle(self.sampleName, self.forwardFastqIn, reverseFastqIn = self.reverseFastqIn, qualityScoreFormat = self.qualityScoreFormat, minWindowQualityAverage = self.minWindowQualityAverage, minLengthAfterTrimming = self.minLengthAfterTrimming, gzippedOutput = self.gzippedOut, clobber = clobber, outputDirectory = outputDirectory)
        self.sickleTrimCommandIndex = jobList.addJob(sickleTrim.sickleCommand)
        sickleTrimJobID = self.submitCommands(mock)
        if self.pairedEndMode:
            self.returnData = {"forward" : runnerSupport.WorkflowReturn(sickleTrim.forwardFastqOut, sickleTrimJobID, sickleTrim.clobber),
                               "reverse" : runnerSupport.WorkflowReturn(sickleTrim.reverseFastqOut, sickleTrimJobID, sickleTrim.clobber),
                               "singleton" : runnerSupport.WorkflowReturn(sickleTrim.singletonFastqOut, sickleTrimJobID, sickleTrim.clobber)}
        else:
            self.returnData = {"forward" : runnerSupport.WorkflowReturn(sickleTrim.forwardFastqOut, sickleTrimJobID, sickleTrim.clobber),
                               "reverse" : runnerSupport.EmptyWorkflowReturn(),
                               "singleton" : runnerSupport.EmptyWorkflowReturn()}
        
    def submitCommands(self, mock):
        import genericRunners
        tempDir = self.tempDir
        sickleTrim = genericRunners.HoffmanJob([self.lastForwardJob.jobID, self.lastReverseJob.jobID], self.sickleTrimCommandIndex, self.sampleName + "Sickle", tempDir, self.emailAddress, "a", 1, mock = mock)
        return sickleTrim.jobID

class OncotatorAnalysis(object):
    
    def __init__(self, jobList, sampleName, somaticVariantPickle = False, oncotatorDB = False, genome = False, emailAddress = False, clobber = False, outputDirectory = "", lastSomaticJob = False, mock = False):
        import houseKeeping
        import gatkRunners
        import programRunners
        import oncotatorRunner
        import runnerSupport
        import os
        self.sampleName = sampleName
        if not somaticVariantPickle and not lastSomaticJob:
            raise RuntimeError("No somatic variant pickle file passed.")
        if somaticVariantPickle:
            self.somaticVariantPickle = somaticVariantPickle
        else:
            self.somaticVariantPickle = lastSomaticJob.data
        if not lastSomaticJob:
            self.lastSomaticJob = runnerSupport.EmptyWorkflowReturn()
        else:
            self.lastSomaticJob = lastSomaticJob
        clobberSet = set([clobber] + [item.clobber for item in [self.lastSomaticJob]])
        if True in clobberSet:
            clobber = True
        elif False  in clobberSet:
            clobber = False
        elif None in clobberSet:
            clobber = None
        self.emailAddress = emailAddress
        self.tempDir = jobList.tempDir
        self.sampleName = sampleName
        oncotatorOutputMaker = programRunners.OncotatorOutputMaker(sampleName, self.somaticVariantPickle, clobber = clobber, outputDirectory = outputDirectory)
        oncotatorRunner = oncotatorRunner.OncotatorWrapper(sampleName, oncotatorOutputMaker.oncotatorMAF, dataBase = oncotatorDB, genome = genome, clobber = oncotatorOutputMaker.clobber, outputDirectory = outputDirectory)
        oncotatorReader = programRunners.OncotatorReader(sampleName, self.somaticVariantPickle, oncotatorRunner.oncotatorOut, oncotatorRunner.clobber, outputDirectory)
        commandList = [oncotatorOutputMaker.oncotatorOutputCommand, oncotatorRunner.oncotatorWrapperCommand, oncotatorReader.oncotatorReaderCommand]
        commandString = " && ".join(commandList)
        self.commandIndex = jobList.addJob(commandString)
        commandJobID = self.submitCommands(mock)
        self.returnData = runnerSupport.WorkflowReturn(oncotatorReader.oncotatorPickle, commandJobID, oncotatorReader.clobber)
        
    def submitCommands(self, mock):
        import genericRunners
        tempDir = self.tempDir
        commandSet = genericRunners.HoffmanJob([self.lastSomaticJob.jobID], self.commandIndex, self.sampleName + "oncotator", tempDir, self.emailAddress, "a", 1, mock = mock)
        return commandSet.jobID

class Athlates(object):
    
    def __init__(self, jobList, sampleName, hlaMolecule, athlatesDB, bamIn = False, emailAddress = False, clobber = False, outputDirectory = "", hlaBAMJob = False, mock = False):
        import houseKeeping
        import gatkRunners
        import programRunners
        import runnerSupport
        import os
        self.sampleName = sampleName
        self.hlaMolecule = hlaMolecule.upper()
        self.athlatesDB = athlatesDB
        if not bamIn and not hlaBAMJob:
            raise RuntimeError("No HLA-aligned BAM file passed.")
        if bamIn:
            self.bamIn = bamIn
        else:
            self.bamIn = hlaBAMJob.data
        if not hlaBAMJob:
            self.hlaBAMJob = runnerSupport.EmptyWorkflowReturn()
        else:
            self.hlaBAMJob = hlaBAMJob
        clobberSet = set([clobber] + [item.clobber for item in [self.hlaBAMJob]])
        if True in clobberSet:
            clobber = True
        elif False  in clobberSet:
            clobber = False
        elif None in clobberSet:
            clobber = None
        self.emailAddress = emailAddress
        self.tempDir = jobList.tempDir
        self.sampleName = sampleName
        bedDir = self.athlatesDB + os.sep + "bed" + os.sep
        self.bedFileOnTarget = bedDir + "hla.%s.bed" %self.hlaMolecule
        self.bedFileOffTarget = bedDir + "hla.non-%s.bed" %self.hlaMolecule
        self.msa = self.athlatesDB + os.sep + "msa" + os.sep + "%s_nuc.txt" %self.hlaMolecule
        self.hlaFasta = self.athlatesDB + os.sep + "ref" + os.sep + "hla.clean.fasta"
        filterAndConvertOnTarget = programRunners.ViewSAMtoBAM(sampleName + ".on", self.bamIn, self.hlaFasta, samInput = False, bamOutput = False, includeHeaders = False, bedFileFilter = self.bedFileOnTarget, clobber = clobber, outputDirectory = outputDirectory)
        filterAndConvertOffTarget = programRunners.ViewSAMtoBAM(sampleName + ".off", self.bamIn, self.hlaFasta, samInput = False, bamOutput = False, includeHeaders = False, bedFileFilter = self.bedFileOffTarget, clobber = filterAndConvertOnTarget.clobber, outputDirectory = outputDirectory)
        sortOnTarget = programRunners.LinuxSort(sampleName, filterAndConvertOnTarget.samOut, self.sampleName + ".on.sorted.sam", ["1,1","3,3"], filterAndConvertOffTarget.clobber, outputDirectory)
        sortOffTarget = programRunners.LinuxSort(sampleName, filterAndConvertOffTarget.samOut, self.sampleName + ".off.sorted.sam", ["1,1","3,3"], sortOnTarget.clobber, outputDirectory)
        compressOnTarget = programRunners.ViewSAMtoBAM(sampleName + ".on.sorted", sortOnTarget.outputFile, self.hlaFasta, samInput = True, bamOutput = True, includeHeaders = False, clobber = sortOffTarget.clobber, outputDirectory = outputDirectory)
        compressOffTarget = programRunners.ViewSAMtoBAM(sampleName + ".off.sorted", sortOffTarget.outputFile, self.hlaFasta, samInput = True, bamOutput = True, includeHeaders = False, clobber = compressOnTarget.clobber, outputDirectory = outputDirectory)
        deleteIntermediates = houseKeeping.Delete([filterAndConvertOnTarget.samOut, filterAndConvertOffTarget.samOut, sortOnTarget.outputFile, sortOffTarget.outputFile])
        athlates = programRunners.Athlates(self.sampleName, compressOnTarget.bamOut, compressOffTarget.bamOut, self.msa, self.hlaMolecule, sortOffTarget.clobber, outputDirectory)
        commandList = [filterAndConvertOnTarget.viewCommand, filterAndConvertOffTarget.viewCommand, sortOnTarget.sortCommand, sortOffTarget.sortCommand, compressOnTarget.viewCommand, compressOffTarget.viewCommand, deleteIntermediates.deleteCommand, athlates.athlatesCommand]
        commandString = " && ".join(commandList)
        self.commandIndex = jobList.addJob(commandString)
        commandJobID = self.submitCommands(mock)
        self.returnData = runnerSupport.WorkflowReturn(athlates.hlaTypingOutput, commandJobID, athlates.clobber)
        
    def submitCommands(self, mock):
        import genericRunners
        tempDir = self.tempDir
        commandSet = genericRunners.HoffmanJob([self.hlaBAMJob.jobID], self.commandIndex, self.sampleName + "HLA" + self.hlaMolecule, tempDir, self.emailAddress, "a", 1, mock = mock)
        return commandSet.jobID
    
class AddHLA(object):
    
    def __init__(self, jobList, sampleName, hlaCallFilesDict, somaticVariantPickle = False, emailAddress = False, clobber = None, outputDir = "", lastSomaticPickleJob = False, hlaJobs = False, mock = False):
        import houseKeeping
        import gatkRunners
        import programRunners
        import runnerSupport
        if not (somaticVariantPickle or lastSomaticPickleJob):
            raise RuntimeError("No somatics job data or somatics file passed.")
        if not hlaCallFilesDict:
            raise RuntimeError("No dictionary of molecule:filepath for HLA call files given.")
        if somaticVariantPickle:
            self.somaticVariantPickle = somaticVariantPickle
        else:
            self.somaticVariantPickle = lastSomaticPickleJob.data
        if not lastSomaticPickleJob:
            lastSomaticPickleJob = runnerSupport.EmptyWorkflowReturn()
        if not hlaJobs:
            hlaJobs = runnerSupport.EmptyWorkflowReturn()
        self.lastSomaticPickleJob = lastSomaticPickleJob
        self.hlaJobs = hlaJobs
        clobberSet = set([clobber] + [item.clobber for item in [self.lastSomaticPickleJob, self.hlaJobs]])
        if True in clobberSet:
            clobber = True
        elif False  in clobberSet:
            clobber = False
        elif None in clobberSet:
            clobber = None
        self.emailAddress = emailAddress
        self.tempDir = jobList.tempDir
        self.sampleName = sampleName
        addHLA = programRunners.AddHLAData(sampleName, self.somaticVariantPickle, hlaCallFilesDict, clobber = clobber, outputDirectory = outputDir)
        self.addHLACommandIndex = jobList.addJob(addHLA.addHLACommand)
        addHLAJobID = self.submitCommands(mock)
        self.returnData = runnerSupport.WorkflowReturn(addHLA.hlaDataOut, addHLAJobID, addHLA.clobber)
        
    def submitCommands(self, mock):
        import genericRunners
        hlaJobIDs = self.hlaJobs.jobID
        if type(hlaJobIDs) == int:
            hlaJobIDs = [hlaJobIDs]
        if self.lastSomaticPickleJob:
            if type(self.lastSomaticPickleJob.jobID) == int:
                hlaJobIDs.append(self.lastSomaticPickleJob.jobID)
            else:
                hlaJobIDs += self.lastSomaticPickleJob.jobID
        tempDir = self.tempDir
        addHLA = genericRunners.HoffmanJob(hlaJobIDs, self.addHLACommandIndex, self.sampleName + "AddHLA", tempDir, self.emailAddress, "a", 1, mock = mock)
        return addHLA.jobID
    
class NetMHCAnalysis(object):
    
    def __init__(self, jobList, sampleName, somaticVariantPickle = False, emailAddress = False, clobber = None, outputDir = "", lastSomaticPickleJob = False, mock = False):
        import houseKeeping
        import variantReaders.netMHCAnalyzer
        import programRunners
        import runnerSupport
        if not (somaticVariantPickle or lastSomaticPickleJob):
            raise RuntimeError("No somatics job data or somatics file passed.")
        if somaticVariantPickle:
            self.somaticVariantPickle = somaticVariantPickle
        else:
            self.somaticVariantPickle = lastSomaticPickleJob.data
        if not lastSomaticPickleJob:
            lastSomaticPickleJob = runnerSupport.EmptyWorkflowReturn()
        self.lastSomaticPickleJob = lastSomaticPickleJob
        clobberSet = set([clobber] + [item.clobber for item in [self.lastSomaticPickleJob]])
        if True in clobberSet:
            clobber = True
        elif False  in clobberSet:
            clobber = False
        elif None in clobberSet:
            clobber = None
        self.emailAddress = emailAddress
        self.tempDir = jobList.tempDir
        self.sampleName = sampleName
        netMHC = variantReaders.netMHCAnalyzer.netMHCWrapper(sampleName, self.somaticVariantPickle, clobber = clobber, outputDirectory = outputDir)
        self.netMHCCommandIndex = jobList.addJob(netMHC.netMHCWrapperCommand)
        netHMCJobID = self.submitCommands(mock)
        self.returnData = runnerSupport.WorkflowReturn(netMHC.peptideTableOut, netHMCJobID, netMHC.clobber)
        
    def submitCommands(self, mock):
        import genericRunners
        tempDir = self.tempDir
        netMHC = genericRunners.HoffmanJob([self.lastSomaticPickleJob.jobID], self.netMHCCommandIndex, self.sampleName + "netMHC", tempDir, self.emailAddress, "a", 1, mock = mock)
        return netMHC.jobID
    
class Capstone(object):
    
    def __init__(self, jobList, sampleName, dependencyList, emailAddress = False, outputDir = "", mock = False):
        import houseKeeping
        self.sampleName = sampleName
        self.dependencyList = dependencyList
        self.tempDir = jobList.tempDir
        self.emailAddress = emailAddress
        capstone = houseKeeping.Capstone(sampleName, outputDir)
        self.capstoneCommandIndex = jobList.addJob(capstone.capstoneCommand)
        self.returnData = self.submitCommands(mock)
        
    def submitCommands(self, mock):
        import genericRunners
        tempDir = self.tempDir
        capstone = genericRunners.HoffmanJob(self.dependencyList, self.capstoneCommandIndex, self.sampleName + "Cap", tempDir, self.emailAddress, "ea", 1, mock = mock)
        return capstone.jobID
