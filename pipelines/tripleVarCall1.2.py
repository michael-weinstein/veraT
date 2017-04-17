#!/usr/bin/env python3

import rootPathSetter
import os
global runnerRoot
runnerRoot = os.sep.join(__file__.split(os.sep)[:-2]) + os.sep

class CheckArgs(object):
    
    def __init__(self):
        import os
        import argparse
        parser = argparse.ArgumentParser()
        parser.add_argument("-m", "--mock", help = "Mock run only, do not submit jobs to scheduler.", action = 'store_true')
        parser.add_argument("-t", "--tumorCriteria", help = "Fields that describe tumor sequence files.", required = True)
        parser.add_argument("-n", "--normalCriteria", help = "Fields that describe normal sequence files.", required = True)
        parser.add_argument("-d", "--sequenceDirectory", help = "Directory containing all sequence files needed (they can be in subdirectories).", required = True)
        parser.add_argument("-r", "--referenceGenomeFasta", help = "Reference genome to use. It should already be BWA and GATK indexed.", default = runnerRoot + "references/UCSC_hg19/ucsc.hg19.fasta")
        parser.add_argument("--hisatIndexData" , help = "Hisat genome index prefix", default = runnerRoot + "references/UCSC_hg19/hisat2/UCSC_hg19")
        parser.add_argument("-i", "--intervals", help = "Intervals file for variant calling.", default = runnerRoot + "references/SeqCap_EZ_Exome_v3_capture_headless.padded100bp.merge.bed")
        parser.add_argument("-j", "--jointGenotypingGVCFDirectory", help = "Directory containing a frozen set of GVCFs for joint genotyping.", default = runnerRoot + "references/HCSNV.hg19")
        parser.add_argument("-e", "--email", help = "Email address to notify for job progress.", default = False)
        parser.add_argument("--dbsnp", help = "dbSNP VCF location", default = runnerRoot + "references/gatkPackage/dbsnp_138.hg19.excluding_sites_after_129.vcf")
        parser.add_argument("--hapmap", help = "Hapmap VCF", default = runnerRoot + "references/gatkPackage/hapmap_3.3.hg19.vcf")
        parser.add_argument("--omni", help = "Omni joint genotyping VCF", default = runnerRoot + "references/gatkPackage/1000G_omni2.5.hg19.sites.vcf")
        parser.add_argument("--snp1000g", help = "1000 genomes SNP VCF", default = runnerRoot + "references/gatkPackage/1000G_phase1.snps.high_confidence.hg19.sites.vcf")
        parser.add_argument("--indel1000g", help = "1000 genomes indel VCF", default = runnerRoot + "references/gatkPackage/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf")        
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
        self.hisatIndexData = rawArgs.hisatIndexData
        
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
    #GATK Haplotype Caller Everything
    haplotypeCallerNormalJob = workFlows.HaplotypeCaller(jobList, normalSampleName, args.refGenomeFasta, intervals = args.intervals, dbSNP = args.dbSNP, emailAddress = args.email, clobber = tumorProcessJobs.clobber, outputDir = outputDir, lastJob = normalProcessJobs, mock = args.mock).returnData
    haplotypeCallerTumorJob = workFlows.HaplotypeCaller(jobList, tumorSampleName, args.refGenomeFasta, intervals = args.intervals, dbSNP = args.dbSNP, emailAddress = args.email, clobber = haplotypeCallerNormalJob.clobber, outputDir = outputDir, lastJob = tumorProcessJobs, mock = args.mock).returnData

    
    #GATK Full joint genotyping on DNA data
    jointGenotypeAndVQSRDNA = workFlows.JointGenotypeAndVQSR(jobList, tumorSampleName, args.refGenomeFasta, args.jointGenotypingGVCFDirectory, args.snpVQSRResources, args.indelVQSRResources, intervals = args.intervals, dbSNP = args.dbSNP, emailAddress = args.email, clobber = haplotypeCallerTumorJob.clobber, outputDir = outputDir, lastJob = [haplotypeCallerNormalJob, haplotypeCallerTumorJob], mock = args.mock).returnData
    
    #Get somatics from GATK
    vcfReaderJob = workFlows.VCFReader(jobList, normalSampleName, tumorSampleName, emailAddress = args.email, clobber = jointGenotypeAndVQSRDNA.clobber, outputDir = outputDir, lastJob = jointGenotypeAndVQSRDNA, mock = args.mock).returnData
    
    #Mutect1
    mutect1SomaticsJob = workFlows.Mutect1ToSomaticsData(jobList, tumorSampleName, args.refGenomeFasta, intervals = args.intervals, dbSNP = args.dbSNP, emailAddress = args.email, clobber = vcfReaderJob.clobber, outputDir = outputDir, lastNormalJob = normalProcessJobs, lastTumorJob = tumorProcessJobs, minDepth = 10, mock = args.mock).returnData

    #VarScan
    varScanSomaticsJob = workFlows.VarScanSomaticsData(jobList, tumorSampleName, args.refGenomeFasta, emailAddress = args.email, clobber = mutect1SomaticsJob.clobber, outputDir = outputDir, lastNormalJob = normalProcessJobs, lastTumorJob = tumorProcessJobs, maxWarnings = 10, minDepth = 10, mock = args.mock).returnData
    
    #VariantCombine
    combineSomaticVariantsJob = workFlows.CombineVarScanHaplotypeCallerMutectSomatics(jobList, tumorSampleName, emailAddress = args.email, clobber = varScanSomaticsJob["indel"].clobber, outputDir = outputDir, varScanCombinedSomaticsJob = varScanSomaticsJob, haplotypeCallerSomaticsJob = vcfReaderJob, mutect1SomaticsJob = mutect1SomaticsJob, mock = args.mock).returnData
  
    fuseTandemVariants = workFlows.combineTandemVariants(jobList, tumorSampleName, emailAddress = args.email, clobber = combineSomaticVariantsJob.clobber, outputDir = outputDir, lastVariantJob = combineSomaticVariantsJob, mock = args.mock).returnData
    
    #Capstone to mark completion
    capstoneDependencies = [fuseTandemVariants.jobID]
    capstoneJobID = workFlows.Capstone(jobList, tumorSampleName, capstoneDependencies, args.email, outputDir, args.mock)
    
    #ending the pipeline
    pipeline.end()

if __name__ == '__main__':
    main()

