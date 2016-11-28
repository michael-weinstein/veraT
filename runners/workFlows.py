#!/usr/bin/env python3

import os

class PipelineStart(object):
    
    def __init__(self, sampleName, tempDir = False):
        if tempDir:
            import os
            if not os.path.isdir(tempDir):
                os.mkdir(tempDir)
            else:
                if os.path.isfile(tempDir + os.sep + "jobs.pkl"):
                    raise RuntimeError("A jobs.pkl already exists in the specified temporary directory.  Please clear the directory before overwriting with the new jobs.")
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
    
    def __init__(self, jobList, sampleName, refGenomeFasta, pe1, pe2 = False, bwaCores = 2, readGroupLibrary = "defaultLibrary", readGroupID = 1, readGroupSampleName = False, readGroupPlatform = "Illumina", barcode = "defaultBarcode", emailAddress = False, clobber = None, outputDir = "", lastJob = False, mock = False):
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
        self.bwaCores = bwaCores
        align = programRunners.BWAmem(sampleName, refGenomeFasta, pe1, pe2, cores = bwaCores, clobber = None, outputDirectory = outputDir)
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
        deleteOriginalBAM = genericRunners.HoffmanJob([addReadGroups.jobID], self.addReadGroupsCommandIndex, self.sampleName + self.library + "AddRG", tempDir, self.emailAddress, "a", cores = 1, mock = mock)
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
        merge = picardRunners.MergeSAMFiles(sampleName, bamFiles, sort_order, validation_stringency, create_index, clobber, outputDir)
        self.mergeCommandIndex = jobList.addJob(merge.mergeSAMCommand)
        self.deleteCommands = []
        for bamFile in bamFiles:
            deletePartialBAM = houseKeeping.Delete(bamFile)
            self.deleteCommands.append(deletePartialBAM.deleteCommand)
        self.deletePartialCommandIndex = jobList.addJob(self.deleteCommands)  #going to have just one node do all the deletions for this job
        deduplicate = picardRunners.Deduplicate(sampleName, merge.bamOut, clobber = merge.clobber, outputDirectory = outputDir)
        self.deduplicateCommandIndex = jobList.addJob(deduplicate.deduplicateCommand)
        deleteOriginalBam = houseKeeping.Delete(merge.bamOut)
        self.deleteOriginalBAMCommandIndex = jobList.addJob(deleteOriginalBam.deleteCommand)
        deduplicateJobID = self.submitCommands(mock)
        self.returnData = runnerSupport.WorkflowReturn(deduplicate.bamOut, deduplicateJobID, deduplicate.clobber)
        
    def submitCommands(self, mock):
        import genericRunners
        tempDir = self.tempDir
        merge = genericRunners.HoffmanJob(self.lastJob.jobID, self.mergeCommandIndex, self.sampleName + "Merge", tempDir, self.emailAddress, "a", 1, mock = mock)
        deletePartialBAM = genericRunners.HoffmanJob(merge.jobID, self.deletePartialCommandIndex, self.sampleName + "DelPart", tempDir, self.emailAddress, "a", cores = 1, mock = mock)
        deduplicate = genericRunners.HoffmanJob(merge.jobID, self.deduplicateCommandIndex, self.sampleName + "Dedup", tempDir, self.emailAddress, "a", cores = 1, mock = mock)
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
        deleteRealignedBAM = genericRunners.HoffmanJob(printReads.jobID, self.deleteRealignedBAMCommandIndex, self.sampleName + "DelRealign", tempDir, self.emailAddress, "a", 1, mock = mock)
        depthOfCoverage = genericRunners.HoffmanJob(printReads.jobID, self.depthOfCoverageCommandIndex, self.sampleName + "depth", tempDir, self.emailAddress, "a", 1, mock = mock)
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
        haplotypeCaller = genericRunners.HoffmanJob(self.lastJob.jobID, self.haplotypeCallerCommandIndex, self.sampleName + "hapCall", tempDir, self.emailAddress, "a", 1, mock = mock)
        jointGenotyper = genericRunners.HoffmanJob(haplotypeCaller.jobID, self.jointGenotypeCommandIndex, self.sampleName + "jointGeno", tempDir, self.emailAddress, "a", 1, mock = mock)
        vqsrSNPAnalysis = genericRunners.HoffmanJob(jointGenotyper.jobID, self.vqsrSNPAnalysisCommandIndex, self.sampleName + "vqsrSNP1", tempDir, self.emailAddress, "a", 1, mock = mock)
        vqsrSNPExecute = genericRunners.HoffmanJob(vqsrSNPAnalysis.jobID, self.vqsrSNPExecuteCommandIndex, self.sampleName + "vqsrSNP2", tempDir, self.emailAddress, "a", 1, mock = mock)
        vqsrIndelAnalysis = genericRunners.HoffmanJob(vqsrSNPExecute.jobID, self.vqsrIndelAnalysisCommandIndex, self.sampleName + "vqsrIndel1", tempDir, self.emailAddress, "a", 1, mock = mock)
        vqsrIndelExecute = genericRunners.HoffmanJob(vqsrIndelAnalysis.jobID, self.vqsrIndelExecuteCommandIndex, self.sampleName + "vqsrIndel2", tempDir, self.emailAddress, "a", 1, mock = mock)
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
    
class VarScanFilter(object):
    
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
        #need to get the script for pulling positions from the data
        #bamReadCount = programRunners.BAMReadCount()
        snpJobID, indelJobID = self.submitCommands(mock)
        self.returnData = {"snp":runnerSupport.WorkflowReturn(varScan.snpOut, snpJobID, varScan.clobber),
                      "indel":runnerSupport.WorkflowReturn(varScan.indelOut, indelJobID, varScan.clobber)}
        
    def submitCommands(self, mock):
        import genericRunners
        tempDir = self.tempDir
        mPileupNormal = genericRunners.HoffmanJob([self.lastNormalJob.jobID], self.mPileupNormalCommandIndex, self.sampleName + "normPile", tempDir, self.emailAddress, "a", 2, mock = mock)
        mPileupTumor = genericRunners.HoffmanJob([self.lastTumorJob.jobID], self.mPileupTumorCommandIndex, self.sampleName + "tumorPile", tempDir, self.emailAddress, "a", 2, mock = mock)
        varScanSomatic = genericRunners.HoffmanJob([mPileupNormal.jobID, mPileupTumor.jobID], self.varScanSomaticCommandIndex, self.sampleName + "varScanSom", tempDir, self.emailAddress, "a", 1, mock = mock)
        varScanSNP = genericRunners.HoffmanJob([varScanSomatic.jobID], self.varScanProcessSNPCommandIndex, self.sampleName + "varScanSNP", tempDir, self.emailAddress, "a", 1, mock = mock)
        varScanIndel = genericRunners.HoffmanJob([varScanSomatic.jobID], self.varScanProcessIndelCommandIndex, self.sampleName  +"varScanIndel", tempDir, self.emailAddress, "a", 1, mock = mock)
        return [varScanIndel.jobID, varScanSNP.jobID]