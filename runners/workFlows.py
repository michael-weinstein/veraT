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
    
    def __init__(self, jobList, sampleName, refGenomeFasta, pe1, pe2 = False, bwaCores = "calculateFromSize", readGroupLibrary = "defaultLibrary", readGroupID = 1, readGroupSampleName = False, readGroupPlatform = "Illumina", barcode = "defaultBarcode", emailAddress = False, clobber = None, outputDir = "", lastJob = False, mock = False):
        import programRunners
        import houseKeeping
        import picardRunners
        import runnerSupport
        if not lastJob:
            lastJob = runnerSupport.EmptyWorkflowReturn()
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
        align = genericRunners.HoffmanJob([], self.bwaCommandIndex, self.sampleName + self.library + "align", tempDir, self.emailAddress, emailStatus, cores = self.bwaCores, mock = mock)
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
            deduplicate = picardRunners.Deduplicate(sampleName, bamFiles[0], clobber = clobber, outputDirectory = outputDir)
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
        self.bqsrAnalysisCommandIndex = jobList.addJob([bqsr.secondPassCommand, bqsr.plotCommand])
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
        findTargets = genericRunners.HoffmanJob(self.lastJob.jobID, self.findTargetsCommandIndex, self.sampleName + "indel1", tempDir, self.emailAddress, "a", 1, mock = mock)
        realigner = genericRunners.HoffmanJob(findTargets.jobID, self.realignerCommandIndex, self.sampleName + "indel2", tempDir, self.emailAddress, "a", 1, mock = mock)
        bqsrFirstPass = genericRunners.HoffmanJob(realigner.jobID, self.bqsrFirstPassCommandIndex, self.sampleName + "bqsr1", tempDir, self.emailAddress, "a", 1, mock = mock)
        printReads = genericRunners.HoffmanJob(bqsrFirstPass.jobID, self.bqsrExecuteCommandIndex, self.sampleName + "bqsr2", tempDir, self.emailAddress, "a", 1, mock = mock)
        indexRecalBAM = genericRunners.HoffmanJob(printReads.jobID, self.indexRecalBAMCommandIndex, self.sampleName + "idxRecal", tempDir, self.emailAddress, "a", 1, mock = mock)
        bqsrAnalysis = genericRunners.HoffmanJob(bqsrFirstPass.jobID, self.bqsrAnalysisCommandIndex, self.sampleName + "bqsrData", tempDir, self.emailAddress, "a", 1, mock = mock)
        deleteRealignedBAM = genericRunners.HoffmanJob([printReads.jobID, bqsrAnalysis.jobID], self.deleteRealignedBAMCommandIndex, self.sampleName + "DelRealign", tempDir, self.emailAddress, "a", 1, mock = mock)
        depthOfCoverage = genericRunners.HoffmanJob([printReads.jobID, indexRecalBAM.jobID], self.depthOfCoverageCommandIndex, self.sampleName + "depth", tempDir, self.emailAddress, "a", 2, mock = mock)
        return indexRecalBAM.jobID
    
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
        haplotypeCaller = gatkRunners.HaplotypeCaller(sampleName, bamFile, refGenomeFasta, intervals, dbSNP = dbSNP, clobber = clobber, outputDirectory = outputDir)
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
        haplotypeCaller = genericRunners.HoffmanJob(self.lastJob.jobID, self.haplotypeCallerCommandIndex, self.sampleName + "hapCall", tempDir, self.emailAddress, "a", 1, memory = 16, mock = mock)
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
        haplotypeCaller = gatkRunners.HaplotypeCaller(sampleName, bamFile, refGenomeFasta, intervals, dbSNP = dbSNP, clobber = clobber, outputDirectory = outputDir)
        self.haplotypeCallerCommandIndex = jobList.addJob(haplotypeCaller.haplotypeCallerCommand)
        haplotypeCallerJobID = self.submitCommands(mock)
        self.returnData = runnerSupport.WorkflowReturn(haplotypeCaller.gvcfOut, haplotypeCallerJobID, haplotypeCaller.clobber)
        
    def submitCommands(self, mock = False):
        import genericRunners
        tempDir = self.tempDir
        haplotypeCaller = genericRunners.HoffmanJob(self.lastJob.jobID, self.haplotypeCallerCommandIndex, self.sampleName + "hapCall", tempDir, self.emailAddress, "a", 1, memory = 16, mock = mock)
        return haplotypeCaller.jobID
    
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
        mutect1 = genericRunners.HoffmanJob([self.lastNormalJob.jobID, self.lastTumorJob.jobID], self.mutect1CommandIndex, self.sampleName + "mutect1", tempDir, self.emailAddress, "a", 1, mock = mock)
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
        normalHaplotypeCaller = HaplotypeCaller(jobList, sampleName, refGenomeFasta, bamFile = normalBAM, intervals = intervals, dbSNP = dbSNP, emailAddress = emailAddress, clobber = clobber, outputDir = outputDir, lastJob = lastNormalJob, mock = mock)
        tumorHaplotypeCaller = HaplotypeCaller(jobList, sampleName, refGenomeFasta, bamFile = tumorBAM, intervals = intervals, dbSNP = dbSNP, emailAddress = emailAddress, clobber = normalHaplotypeCaller.returnData.clobber, outputDir = outputDir, lastJob = lastTumorJob, mock = mock)
        jointGenotyping = JointGenotypeAndVQSR(jobList, sampleName, refGenomeFasta, jointGenotypingFrozenGVCFDirectory, snpResources, indelResources, gvcfs = [normalHaplotypeCaller.returnData.data, tumorHaplotypeCaller.returnData.data], intervals = intervals, dbSNP = dbSNP, emailAddress = emailAddress, clobber = tumorHaplotypeCaller.returnData.clobber, outputDir = outputDir, lastJob = [normalHaplotypeCaller.returnData, tumorHaplotypeCaller.returnData], mock = mock)
        vcfReader = programRunners.VCFReader(sampleName, jointGenotyping.returnData.data, tumorBaseName, normalBaseName, maxPValue, minDepth, jointGenotyping.returnData.clobber)
        self.vcfReaderJobIndex = jobList.addJob(vcfReader.vcfReaderCommand)
        vcfReaderJobID = self.submitCommands(mock, jointGenotyping)
        self.returnData = runnerSupport.WorkflowReturn(vcfReader.acceptedOut, vcfReaderJobID, vcfReader.clobber)
        
    def submitCommands(self, mock, jointGenotyping):
        import genericRunners
        tempDir = self.tempDir
        vcfReader = genericRunners.HoffmanJob([jointGenotyping.returnData.jobID], self.vcfReaderJobIndex, self.sampleName + "VCFRead", tempDir, self.emailAddress, "a", 1, mock = mock)
        return vcfReader.jobID
    
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
        mPileupNormal = programRunners.MPileup(sampleName + "norm", normalBAM, refGenomeFasta, bgzip = True, clobber = clobber, outputDirectory = outputDir)
        self.mPileupNormalCommandIndex = jobList.addJob(mPileupNormal.mPileupCommand)
        mPileupTumor = programRunners.MPileup(sampleName + "tumor", tumorBAM, refGenomeFasta, bgzip = True, clobber = mPileupNormal.clobber, outputDirectory = outputDir)
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
        mPileupNormal = genericRunners.HoffmanJob([self.lastNormalJob.jobID], self.mPileupNormalCommandIndex, self.sampleName + "normPile", tempDir, self.emailAddress, "a", 2, mock = mock)
        mPileupTumor = genericRunners.HoffmanJob([self.lastTumorJob.jobID], self.mPileupTumorCommandIndex, self.sampleName + "tumorPile", tempDir, self.emailAddress, "a", 2, mock = mock)
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
        snpVarScanReader = programRunners.VarScanReader(self, sampleName, snpFPFilter.passOut, minDepth = minDepth, clobber = clobber, outputDirectory = outputDir)
        self.snpVarScanReaderCommandIndex = jobList.addJob(snpVarScanReader.varScanReaderCommand)
        indelVarScanReader = programRunners.VarScanReader(self, sampleName, indelFPFilter.passOut, minDepth = minDepth, clobber = clobber, outputDirectory = outputDir)
        self.indelVarScanReaderCommandIndex = jobList.addJob(snpVarScanReader.varScanReaderCommand)
        snpJobID, indelJobID = self.submitCommands(mock, varScan)
        self.returnData = {"snp":runnerSupport.WorkflowReturn(snpVarScanReader.acceptedOut, snpJobID, snpVarScanReader.clobber),
                           "indel":runnerSupport.WorkflowReturn(indelVarScanReader.acceptedOut, snpJobID, indelVarScanReader.clobber)}
        
    def submitCommands(self, mock, varScan):
        import genericRunners
        tempDir = self.tempDir
        snpPositionPuller = genericRunners.HoffmanJob([varScan.returnData["snp"].jobID], self.snpPositionPullerCommandIndex, self.sampleName + "snpPos", tempDir, self.emailAddress, "a", 1, mock = mock)
        indelPositionPuller = genericRunners.HoffmanJob([varScan.returnData["indel"].jobID], self.indelPositionPullerCommandIndex, self.sampleName + "indelPos", tempDir, self.emailAddress, "a", 1, mock = mock)
        snpBAMReadCount = genericRunners.HoffmanJob([snpPositionPuller.jobID], self.snpBAMReadCountCommandIndex, self.sampleName + "snpBAMRead", tempDir, self.emailAddress, "a", 1, mock = mock)
        indelBAMReadCount = genericRunners.HoffmanJob([indelPositionPuller.jobID], self.indelBAMReadCountCommandIndex, self.sampleName + "indelBAMRead", tempDir, self.emailAddress, "a", 1, mock = mock)
        snpFPFilter = genericRunners.HoffmanJob([snpBAMReadCount.jobID], self.snpFPFilterCommandIndex, self.sampleName + "snpFPF", tempDir, self.emailAddress, "a", 1, mock = mock)
        indelFPFilter = genericRunners.HoffmanJob([indelBAMReadCount.jobID], self.indelFPFilterCommandIndex, self.sampleName + "indelFPF", tempDir, self.emailAddress, "a", 1, mock = mock)
        snpVarScanReader = genericRunners.HoffmanJob([snpFPFilter.jobID], self.snpFPFilterCommandIndex, self.sampleName + "snpVSR", tempDir, self.emailAddress, "a", 1, mock = mock)
        indelVarScanReader = genericRunners.HoffmanJob([indelFPFilter.jobID], self.indelVarScanReaderCommandIndex, self.sampleName + "indelVSR", tempDir, self.emailAddress, "a", 1, mock = mock)
        return [snpVarScanReader.jobID, indelVarScanReader.jobID]
    
class CombineVarScanHaplotypeCallerMutectSomatics(object):
    
    def __init__(self, jobList, sampleName, varScanSomaticIndels = False, varScanSomaticSNPs = False, haplotypeCallerSomatics = False, mutect1Somatics = False, emailAddress = False, clobber = None, outputDir = "", varScanSNPSomaticsJob = False, varScanIndelSomaticsJob = False, varScanCombinedSomaticsJob = False, haplotypeCallerSomaticsJob = False, mutect1SomaticsJob = False, mock = False):
        import houseKeeping
        import gatkRunners
        import programRunners
        import runnerSupport
        if not (varScanSomaticIndels and varScanSomaticSNPs) or (varScanIndelSomaticsJob and varScanSNPSomaticsJob) or varScanCombinedSomaticsJob:
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
        if not tumorBAM and not lastTumorJob:
            raise RuntimeError("No tumor bam for processing specified in last job data or arguments.")
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
        self.returnData = runnerSupport.WorkflowReturn(variantCombine.acceptedSomaticVariants, vcfReaderJobID, variantCombine.clobber)
        
    def submitCommands(self, mock):
        import genericRunners
        tempDir = self.tempDir
        variantCombine = genericRunners.HoffmanJob([self.varScanIndelSomaticsJob.jobID, self.varScanSNPSomaticsJob.jobID, self.haplotypeCallerSomaticsJob.jobID, self.mutect1SomaticsJob.jobID], self.variantCombineCommandIndex, self.sampleName + "Combine", tempDir, self.emailAddress, "a", 1, mock = mock)
        return variantCombine.jobID
    
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
