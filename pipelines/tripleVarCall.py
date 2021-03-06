#!/usr/bin/env python3

# class CheckArgsOld(object):
#     
#     def __init__(self):
#         self.mock = False
#         self.tumorCriteria = ["JMS", "scalp"]
#         self.normalCriteria = ["JMS", "normal"]
#         self.sequenceDirectory = "sampleData"
#         self.refGenomeFasta = "references/UCSC_hg19/ucsc.hg19.fasta"
#         self.intervals = "references/SeqCap_EZ_Exome_v3_capture_headless.padded100bp.merge.bed"
#         self.jointGenotypingGVCFDirectory = "references/HCSNV.hg19"
#         self.email = "mweinste@ucla.edu"
#         self.dbSNP = "references/gatkPackage/dbsnp_138.hg19.excluding_sites_after_129.vcf"
#         self.hapmap = "references/gatkPackage/hapmap_3.3.hg19.vcf"
#         self.omni = "references/gatkPackage/1000G_omni2.5.hg19.sites.vcf"
#         self.snp1000g = "references/gatkPackage/1000G_phase1.snps.high_confidence.hg19.sites.vcf"
#         self.indel1000g = "references/gatkPackage/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"
#         self.makeRecalRefs()
#         
#     def makeRecalRefs(self):
#         from runners.gatkRunners import VQSRResource
#         #making SNP resource list
#         self.snpVQSRResources = []
#         self.snpVQSRResources.append(VQSRResource(self.dbSNP, "dbsnp", known = True, training = False, truth = False, prior = 2.0))
#         self.snpVQSRResources.append(VQSRResource(self.omni, "omni", False, True, True, 12.0))
#         self.snpVQSRResources.append(VQSRResource(self.snp1000g, "1000G", False, True, False, 10.0))
#         self.snpVQSRResources.append(VQSRResource(self.hapmap, "hapmap", False, True, True, 15.0))
# 
#         #making indel resource list
#         self.indelVQSRResources = []
#         self.indelVQSRResources.append(VQSRResource(self.indel1000g, "mills", False, True, True, 12.0))
#         self.indelVQSRResources.append(VQSRResource(self.dbSNP, "dbsnp", True, False, False, 2.0))

class CheckArgs(object):
    
    def __init__(self):
        import os
        import argparse
        parser = argparse.ArgumentParser()
        parser.add_argument("-m", "--mock", help = "Mock run only, do not submit jobs to scheduler.", action = 'store_true')
        parser.add_argument("-t", "--tumorCriteria", help = "Fields that describe tumor sequence files.", required = True)
        parser.add_argument("-n", "--normalCriteria", help = "Fields that describe normal sequence files.", required = True)
        parser.add_argument("-d", "--sequenceDirectory", help = "Directory containing all sequence files needed (they can be in subdirectories).", required = True)
        parser.add_argument("-r", "--referenceGenomeFasta", help = "Reference genome to use. It should already be BWA and GATK indexed.", default = "references/UCSC_hg19/ucsc.hg19.fasta")
        parser.add_argument("-i", "--intervals", help = "Intervals file for variant calling.", default = "references/SeqCap_EZ_Exome_v3_capture_headless.padded100bp.merge.bed")
        parser.add_argument("-j", "--jointGenotypingGVCFDirectory", help = "Directory containing a frozen set of GVCFs for joint genotyping.", default = "references/HCSNV.hg19")
        parser.add_argument("-e", "--email", help = "Email address to notify for job progress.", default = False)
        parser.add_argument("--dbsnp", help = "dbSNP VCF location", default = "references/gatkPackage/dbsnp_138.hg19.excluding_sites_after_129.vcf")
        parser.add_argument("--hapmap", help = "Hapmap VCF", default = "references/gatkPackage/hapmap_3.3.hg19.vcf")
        parser.add_argument("--omni", help = "Omni joint genotyping VCF", default = "references/gatkPackage/1000G_omni2.5.hg19.sites.vcf")
        parser.add_argument("--snp1000g", help = "1000 genomes SNP VCF", default = "references/gatkPackage/1000G_phase1.snps.high_confidence.hg19.sites.vcf")
        parser.add_argument("--indel1000g", help = "1000 genomes indel VCF", default = "references/gatkPackage/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf")        
        parser.add_argument("-o", "--outputDirectory", help = "Directory for data output", required = True)
        rawArgs = parser.parse_args()
        self.mock = rawArgs.mock
        self.tumorCriteria = rawArgs.tumorCriteria.split(",")
        self.normalCriteria = rawArgs.normalCriteria.split(",")
        sequenceDirectory = rawArgs.sequenceDirectory
        if not os.path.isdir(sequenceDirectory):
            raise RuntimeError("Sequence Directory %s not found." %(sequenceDirectory))
        self.sequenceDirectory = sequenceDirectory
        outputDir = rawArgs.outputDirectory
        if not os.path.isdir(outputDir):
            os.mkdir(outputDir)
        self.outputDir = outputDir
        self.refGenomeFasta = rawArgs.referenceGenomeFasta
        self.intervals = rawArgs.intervals
        self.jointGenotypingGVCFDirectory = rawArgs.jointGenotypingGVCFDirectory
        self.email = rawArgs.email
        self.dbSNP = rawArgs.dbsnp
        self.hapmap = rawArgs.hapmap
        self.omni = rawArgs.omni
        self.snp1000g = rawArgs.snp1000g
        self.indel1000g = rawArgs.indel1000g
        self.makeRecalRefs()
        
    def makeRecalRefs(self):
        from runners.gatkRunners import VQSRResource
        #making SNP resource list
        self.snpVQSRResources = []
        self.snpVQSRResources.append(VQSRResource(self.dbSNP, "dbsnp", known = True, training = False, truth = False, prior = 2.0))
        self.snpVQSRResources.append(VQSRResource(self.omni, "omni", False, True, True, 12.0))
        self.snpVQSRResources.append(VQSRResource(self.snp1000g, "1000G", False, True, False, 10.0))
        self.snpVQSRResources.append(VQSRResource(self.hapmap, "hapmap", False, True, True, 15.0))

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
    tumorSampleName = list(tumorSeqTree.keys())[0]
    outputDir = args.outputDir
    normalSampleName = list(normalSeqTree.keys())[0]
    patientID = args.tumorCriteria[0]  #could also get this from normalCriteria
    pipeline = workFlows.PipelineStart(tumorSampleName, useDir = outputDir)
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
                normalAlignJobs = workFlows.AlignMemAddReadGroups(jobList, normalSampleName + str(lane), args.refGenomeFasta, currentLane[1].filePath, currentLane[2].filePath, readGroupSampleName = normalSampleName, readGroupID = lane, clobber = clobber, emailAddress = args.email, outputDir = outputDir, mock = args.mock).returnData
            else:  #if the job list has been initiated here (so it's not false), add the jobs to it
                normalAlignJobs.addJob(workFlows.AlignMemAddReadGroups(jobList, normalSampleName + str(lane), args.refGenomeFasta, currentLane[1].filePath, currentLane[2].filePath, readGroupSampleName = normalSampleName, readGroupID = lane, clobber = normalAlignJobs.clobber, emailAddress = args.email, outputDir = outputDir, mock = args.mock).returnData)
    
    #iteratively set up alignment jobs for tumor sample (iterate over sample, then lane, then take pairs of paired ends if available)
    tumorAlignJobs = False
    for sample in keys(tumorSeqTree): #should be only one here
        for lane in keys(tumorSeqTree[sample]):
            currentLane = tumorSeqTree[sample][lane]
            if not 2 in currentLane:
                currentLane[2] = FastqFile(False)
            sampleFile = currentLane[1]
            if not tumorAlignJobs:
                tumorAlignJobs = workFlows.AlignMemAddReadGroups(jobList, tumorSampleName + str(lane), args.refGenomeFasta, currentLane[1].filePath, currentLane[2].filePath, readGroupSampleName = tumorSampleName, readGroupID = lane, clobber = normalAlignJobs.clobber, emailAddress = args.email, outputDir = outputDir, mock = args.mock).returnData
            else:
                tumorAlignJobs.addJob(workFlows.AlignMemAddReadGroups(jobList, tumorSampleName + str(lane), args.refGenomeFasta, currentLane[1].filePath, currentLane[2].filePath, readGroupSampleName = tumorSampleName, readGroupID = lane, clobber = tumorAlignJobs.clobber, emailAddress = args.email, outputDir = outputDir, mock = args.mock).returnData)

    #merge and deduplicate normal and tumor samples
    normalDeduplicateJobs = workFlows.MergeBAMFilesAndDeduplicate(jobList, normalSampleName, clobber = tumorAlignJobs.clobber, emailAddress = args.email, outputDir = outputDir, lastJob = normalAlignJobs, mock = args.mock).returnData
    tumorDeduplicateJobs = workFlows.MergeBAMFilesAndDeduplicate(jobList, tumorSampleName, clobber = normalDeduplicateJobs.clobber, emailAddress = args.email, outputDir = outputDir, lastJob = tumorAlignJobs, mock = args.mock).returnData
    
    #GATK process BAMs
    normalProcessJobs = workFlows.GATKProcessBAM(jobList, normalSampleName, args.refGenomeFasta, dbSNP = args.dbSNP, intervals = args.intervals, clobber = tumorDeduplicateJobs.clobber, emailAddress = args.email, outputDir = outputDir, lastJob = normalDeduplicateJobs, mock = args.mock).returnData
    tumorProcessJobs = workFlows.GATKProcessBAM(jobList, tumorSampleName, args.refGenomeFasta, dbSNP = args.dbSNP, intervals = args.intervals, clobber = normalProcessJobs.clobber, emailAddress = args.email, outputDir = outputDir, lastJob = tumorDeduplicateJobs, mock = args.mock).returnData
    
    #Run variant calling
    #GATK Haplotype Caller
    normalHCJobs = workFlows.HaplotypeCallerAndVQSR(jobList, normalSampleName, args.refGenomeFasta, args.jointGenotypingGVCFDirectory, args.snpVQSRResources, args.indelVQSRResources, intervals = args.intervals, dbSNP = args.dbSNP, emailAddress = args.email, clobber = tumorProcessJobs.clobber, outputDir = outputDir, lastJob = normalProcessJobs, mock = args.mock).returnData
    tumorHCJobs = workFlows.HaplotypeCallerAndVQSR(jobList, tumorSampleName, args.refGenomeFasta, args.jointGenotypingGVCFDirectory, args.snpVQSRResources, args.indelVQSRResources, intervals = args.intervals, dbSNP = args.dbSNP, emailAddress = args.email, clobber = tumorProcessJobs.clobber, outputDir = outputDir, lastJob = tumorProcessJobs, mock = args.mock).returnData
    #Mutect1
    mutect1Jobs = workFlows.Mutect1(jobList, tumorSampleName, args.refGenomeFasta, intervals = args.intervals, dbSNP = args.dbSNP, emailAddress = args.email, clobber = tumorHCJobs.clobber, outputDir = outputDir, lastNormalJob = normalProcessJobs, lastTumorJob = tumorProcessJobs, mock = args.mock).returnData
    #VarScan
    varScanJobs = workFlows.VarScanFilter(jobList, tumorSampleName, args.refGenomeFasta, emailAddress = args.email, clobber = mutect1Jobs.clobber, outputDir = outputDir, lastNormalJob = normalProcessJobs, lastTumorJob = tumorProcessJobs, mock = args.mock).returnData

    #Capstone to mark completion
    capstoneDependencies = [varScanJobs["snp"].jobID, varScanJobs["indel"].jobID, mutect1Jobs.jobID, normalHCJobs.jobID, tumorHCJobs.jobID]
    capstoneJobID = workFlows.Capstone(jobList, tumorSampleName, capstoneDependencies, args.email, outputDir, args.mock)
    
    #ending the pipeline
    pipeline.end()

if __name__ == '__main__':
    main()

