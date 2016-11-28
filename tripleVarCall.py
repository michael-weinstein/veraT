#!/usr/bin/env python3

import shutil
import os
if os.path.isdir("JMS"):
    shutil.rmtree("JMS")

class CheckArgs(object):
    
    def __init__(self):
        self.mock = True
        self.tumorCriteria = ["JMS", "scalp"]
        self.normalCriteria = ["JMS", "normal"]
        self.sequenceDirectory = "sampleData"
        self.refGenomeFasta = "references/UCSC_hg19/ucsc.hg19.fasta"
        self.intervals = "references/SeqCap_EZ_Exome_v3_capture_headless.padded100bp.merge.bed"
        self.jointGenotypingGVCFDirectory = "references/HCSNV.hg19"
        self.outputDir = self.tumorCriteria[0]
        import os
        if not os.path.isdir(self.outputDir):
            os.mkdir(self.outputDir)
        self.email = "mweinste@ucla.edu"
        self.dbSNP = "references/gatkPackage/dbsnp_138.hg19.excluding_sites_after_129.vcf"
        self.hapmap = None  #hapmap is no longer included with the GATK package download.  Might need to copy from your system to match versions
        self.omni = "references/gatkPackage/1000G_omni2.5.hg19.sites.vcf"
        self.snp1000g = "references/gatkPackage/1000G_phase1.snps.high_confidence.hg19.sites.vcf"
        self.indel1000g = "references/gatkPackage/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"
        self.makeRecalRefs()
        
    def makeRecalRefs(self):
        from runners.gatkRunners import VQSRResource
        #making SNP resource list
        self.snpVQSRResources = []
        self.snpVQSRResources.append(VQSRResource(self.dbSNP, "dbsnp", known = True, training = False, truth = False, prior = 2.0))
        self.snpVQSRResources.append(VQSRResource(self.omni, "omni", False, True, True, 12.0))
        self.snpVQSRResources.append(VQSRResource(self.snp1000g, "1000G", False, True, False, 10.0))
        #making indel resource list
        self.indelVQSRResources = []
        self.indelVQSRResources.append(VQSRResource(self.indel1000g, "mills", False, True, True, 12.0))
        self.indelVQSRResources.append(VQSRResource(self.dbSNP, "dbsnp", True, False, False, 2.0))
 
def getFileMatrix(directory, criteria):
    from runners import fastqDirectoryParser
    parser = fastqDirectoryParser.FastqDirectory(directory)
    return parser.getAllFilesForSample(criteria)

def keys(dictionary):
    assert type(dictionary) == dict
    return list(dictionary.keys())

def main():
    global args
    args = CheckArgs()
    clobber = None
    from runners import workFlows
    from runners.fastqDirectoryParser import FastqFile
    normalSeqTree = getFileMatrix(args.sequenceDirectory, args.normalCriteria)
    tumorSeqTree = getFileMatrix(args.sequenceDirectory, args.tumorCriteria)
    sampleName = args.tumorCriteria[0]  #could also get this from normalCriteria
    pipeline = workFlows.PipelineStart(sampleName)
    jobList = pipeline.getJobList()  #initialize the pipeline and job list
    
    #iteratively set up alignment jobs for nomral tissue
    normalAlignJobs = False
    for sample in keys(normalSeqTree): #should be only one here
        for lane in keys(normalSeqTree[sample]):
            currentLane = normalSeqTree[sample][lane]
            if not 2 in currentLane: #which is to say that we did not have a paired end read here
                currentLane[2] = FastqFile(False)
            sampleFile = currentLane[1] #should be a fastqFile object, could use either one, but not guaranteed a second paired end
            if not normalAlignJobs: #if we have not initiated the job list for this yet, do that
                normalAlignJobs = workFlows.AlignMemAddReadGroups(jobList, sampleName + "norm" + str(lane), args.refGenomeFasta, currentLane[1].filePath, currentLane[2].filePath, readGroupSampleName = sampleName + "norm", readGroupID = lane, clobber = clobber, emailAddress = args.email, outputDir = args.outputDir, mock = args.mock).returnData
            else:  #if the job list has been initiated here (so it's not false), add the jobs to it
                normalAlignJobs.addJob(workFlows.AlignMemAddReadGroups(jobList, sampleName + "norm" + str(lane), args.refGenomeFasta, currentLane[1].filePath, currentLane[2].filePath, readGroupSampleName = sampleName + "norm", readGroupID = lane, clobber = normalAlignJobs.clobber, emailAddress = args.email, outputDir = args.outputDir, mock = args.mock).returnData)
    
    #iteratively set up alignment jobs for tumor sample (iterate over sample, then lane, then take pairs of paired ends if available)
    tumorAlignJobs = False
    for sample in keys(tumorSeqTree): #should be only one here
        for lane in keys(tumorSeqTree[sample]):
            currentLane = tumorSeqTree[sample][lane]
            if not 2 in currentLane:
                currentLane[2] = FastqFile(False)
            sampleFile = currentLane[1]
            if not tumorAlignJobs:
                tumorAlignJobs = workFlows.AlignMemAddReadGroups(jobList, sampleName + "tumor" + str(lane), args.refGenomeFasta, currentLane[1].filePath, currentLane[2].filePath, readGroupSampleName = sampleName + "tumor", readGroupID = lane, clobber = normalAlignJobs.clobber, emailAddress = args.email, outputDir = args.outputDir, mock = args.mock).returnData
            else:
                tumorAlignJobs.addJob(workFlows.AlignMemAddReadGroups(jobList, sampleName + "tumor" + str(lane), args.refGenomeFasta, currentLane[1].filePath, currentLane[2].filePath, readGroupSampleName = sampleName + "tumor", readGroupID = lane, clobber = tumorAlignJobs.clobber, emailAddress = args.email, outputDir = args.outputDir, mock = args.mock).returnData)
                
    #merge and deduplicate normal and tumor samples
    normalDeduplicateJobs = workFlows.MergeBAMFilesAndDeduplicate(jobList, sampleName + "norm", clobber = tumorAlignJobs.clobber, emailAddress = args.email, outputDir = args.outputDir, lastJob = normalAlignJobs, mock = args.mock).returnData
    tumorDeduplicateJobs = workFlows.MergeBAMFilesAndDeduplicate(jobList, sampleName + "tumor", clobber = normalDeduplicateJobs.clobber, emailAddress = args.email, outputDir = args.outputDir, lastJob = tumorAlignJobs, mock = args.mock).returnData
    
    #GATK process BAMs
    normalProcessJobs = workFlows.GATKProcessBAM(jobList, sampleName + "norm", args.refGenomeFasta, dbSNP = args.dbSNP, clobber = tumorDeduplicateJobs.clobber, emailAddress = args.email, outputDir = args.outputDir, lastJob = normalDeduplicateJobs, mock = args.mock).returnData
    tumorProcessJobs = workFlows.GATKProcessBAM(jobList, sampleName, args.refGenomeFasta, dbSNP = args.dbSNP, clobber = normalProcessJobs.clobber, emailAddress = args.email, outputDir = args.outputDir, lastJob = tumorDeduplicateJobs, mock = args.mock).returnData
    
    #Run variant calling
    #GATK Haplotype Caller
    normalHCJobs = workFlows.HaplotypeCallerAndVQSR(jobList, sampleName + "norm", args.refGenomeFasta, args.jointGenotypingGVCFDirectory, args.snpVQSRResources, args.indelVQSRResources, intervals = args.intervals, dbSNP = args.dbSNP, emailAddress = args.email, clobber = tumorProcessJobs.clobber, outputDir = args.outputDir, lastJob = normalProcessJobs, mock = args.mock).returnData
    tumorHCJobs = workFlows.HaplotypeCallerAndVQSR(jobList, sampleName + "tumor", args.refGenomeFasta, args.jointGenotypingGVCFDirectory, args.snpVQSRResources, args.indelVQSRResources, intervals = args.intervals, dbSNP = args.dbSNP, emailAddress = args.email, clobber = tumorProcessJobs.clobber, outputDir = args.outputDir, lastJob = tumorProcessJobs, mock = args.mock).returnData
    #Mutect1
    mutect1Jobs = workFlows.Mutect1(jobList, sampleName, args.refGenomeFasta, intervals = args.intervals, dbSNP = args.dbSNP, emailAddress = args.email, clobber = tumorHCJobs.clobber, outputDir = args.outputDir, lastNormalJob = normalProcessJobs, lastTumorJob = tumorProcessJobs, mock = args.mock).returnData
    #VarScan
    varScanJobs = workFlows.VarScanFilter(jobList, sampleName, args.refGenomeFasta, emailAddress = args.email, clobber = mutect1Jobs.clobber, outputDir = args.outputDir, lastNormalJob = normalProcessJobs, lastTumorJob = tumorProcessJobs, mock = args.mock).returnData

    #ending the pipeline
    pipeline.end()

if __name__ == '__main__':
    main()

