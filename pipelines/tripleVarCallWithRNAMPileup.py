#!/usr/bin/env python3

import rootPathSetter
import os
global runnerRoot
runnerRoot = os.sep.join(os.path.abspath(__file__).split(os.sep)[:-2]) + os.sep

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
        parser.add_argument("--athlatesDB", help = "HLA reference genome to use (probably part of athlates package).", default = runnerRoot + "references/athlatesDB/")
        parser.add_argument("--hisatIndexData" , help = "Hisat genome index prefix", default = runnerRoot + "references/UCSC_hg19/hisat2/UCSC_hg19")
        parser.add_argument("-a", "--rnaCriteria", help = "Fields that describe RNA sequence data", required = True)
        parser.add_argument("-i", "--intervals", help = "Intervals file for variant calling.", default = runnerRoot + "references/SeqCap_EZ_Exome_v3_capture_headless.padded100bp.merge.bed")
        parser.add_argument("-j", "--jointGenotypingGVCFDirectory", help = "Directory containing a frozen set of GVCFs for joint genotyping.", default = runnerRoot + "references/HCSNV.hg19")
        parser.add_argument("-e", "--email", help = "Email address to notify for job progress.", default = False)
        parser.add_argument("--dbsnp", help = "dbSNP VCF location", default = runnerRoot + "references/gatkPackage/dbsnp_138.hg19.excluding_sites_after_129.vcf")
        parser.add_argument("--hapmap", help = "Hapmap VCF", default = runnerRoot + "references/gatkPackage/hapmap_3.3.hg19.vcf")
        parser.add_argument("--omni", help = "Omni joint genotyping VCF", default = runnerRoot + "references/gatkPackage/1000G_omni2.5.hg19.sites.vcf")
        parser.add_argument("--snp1000g", help = "1000 genomes SNP VCF", default = runnerRoot + "references/gatkPackage/1000G_phase1.snps.high_confidence.hg19.sites.vcf")
        parser.add_argument("--indel1000g", help = "1000 genomes indel VCF", default = runnerRoot + "references/gatkPackage/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf")        
        parser.add_argument("-o", "--outputDirectory", help = "Directory for data output", required = True)
        parser.add_argument("--peptideDictionaryPickle", help = "Pickled dictionary of protein sequence and ID data.", default = runnerRoot + "/references/peptides/Homo_sapiens.GRCh37.75.pep.all.fa.pkl")
        rawArgs = parser.parse_args()
        self.mock = rawArgs.mock
        self.tumorCriteria = rawArgs.tumorCriteria.split(",")
        self.normalCriteria = rawArgs.normalCriteria.split(",")
        self.rnaCriteria = rawArgs.rnaCriteria.split(",")
        sequenceDirectory = rawArgs.sequenceDirectory
        if not os.path.isdir(sequenceDirectory):
            raise RuntimeError("Sequence Directory %s not found." %(sequenceDirectory))
        self.sequenceDirectory = sequenceDirectory
        outputDir = rawArgs.outputDirectory
        if not os.path.isdir(outputDir):
            os.mkdir(outputDir)
        self.outputDir = outputDir
        self.refGenomeFasta = rawArgs.referenceGenomeFasta
        self.athlatesDB = rawArgs.athlatesDB
        peptideDictionaryPickle = rawArgs.peptideDictionaryPickle
        if not os.path.isfile(peptideDictionaryPickle):
            raise FileNotFoundError("Unable to find peptide dictionary pickle at %s" %peptideDictionaryPickle)
        self.peptideDictionaryPickle = peptideDictionaryPickle
        self.hlaFasta = self.athlatesDB + "/ref/hla.clean.fasta"
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
    rnaSeqTree = getFileMatrix(args.sequenceDirectory, args.rnaCriteria)
    tumorSampleName = list(tumorSeqTree.keys())[0]
    outputDir = args.outputDir
    normalSampleName = list(normalSeqTree.keys())[0]
    patientID = args.tumorCriteria[0]  #could also get this from normalCriteria
    pipeline = workFlows.PipelineStart(tumorSampleName, useDir = outputDir)
    jobList = pipeline.getJobList()  #initialize the pipeline and job list

#SOMATIC VARIANT FINDING SECTION
    #iteratively set up alignment jobs for normal tissue
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
    
    #Combine outputs from the three variant callers
    combineSomaticVariantsJob = workFlows.CombineVarScanHaplotypeCallerMutectSomatics(jobList, tumorSampleName, emailAddress = args.email, clobber = varScanSomaticsJob["indel"].clobber, outputDir = outputDir, varScanCombinedSomaticsJob = varScanSomaticsJob, haplotypeCallerSomaticsJob = vcfReaderJob, mutect1SomaticsJob = mutect1SomaticsJob, mock = args.mock).returnData

#RNA ANALYSIS SECTION
        # starts with iteratively set up alignment jobs for tumor RNA sample (iterate over sample, then lane, then take pairs of paired ends if available)
    rnaAlignJobs = False
    #print ("Keys: %s" %(keys(rnaSeqTree)))
    for sample in keys(rnaSeqTree): #should be only one here
        for lane in keys(rnaSeqTree[sample]):
            currentLane = rnaSeqTree[sample][lane]
            if not 2 in currentLane:
                currentLane[2] = FastqFile(False)
            sampleFile = currentLane[1]
            if not rnaAlignJobs:
                #print("Hit the first")
                rnaAlignJobs = workFlows.AlignHisat2AddReadGroups(jobList, tumorSampleName + "RNA" + str(lane), args.refGenomeFasta, args.hisatIndexData, currentLane[1].filePath, pe2 = currentLane[2].filePath, readGroupID = lane, readGroupSampleName = tumorSampleName + "RNA", emailAddress = False, clobber = combineSomaticVariantsJob.clobber, outputDir = outputDir, mock = args.mock).returnData
            else:
                #print("Hit the second")
                rnaAlignJobs.addJob(workFlows.AlignHisat2AddReadGroups(jobList, tumorSampleName + "RNA" + str(lane), args.refGenomeFasta, args.hisatIndexData, currentLane[1].filePath, pe2 = currentLane[2].filePath, readGroupID = lane, readGroupSampleName = tumorSampleName + "RNA", emailAddress = False, clobber = rnaAlignJobs.clobber, outputDir = outputDir, mock = args.mock).returnData)
    #print(rnaAlignJobs)
    if not rnaAlignJobs:
        raise RuntimeError("Got back no RNA alignment jobs. Please check your definition of RNA sequence data files.")
    rnaDeduplicateJobs = workFlows.MergeBAMFilesAndDeduplicate(jobList, tumorSampleName + "RNA", clobber = rnaAlignJobs.clobber, emailAddress = args.email, outputDir = outputDir, lastJob = rnaAlignJobs, mock = args.mock).returnData
    rnaSplitAndTrim = workFlows.SplitAndTrimRNAData(jobList, tumorSampleName + "RNA", args.refGenomeFasta, emailAddress = args.email, clobber = rnaDeduplicateJobs.clobber, outputDir = outputDir, lastJob = rnaDeduplicateJobs, mock = args.mock).returnData
    rnaProcessJobs = workFlows.GATKProcessBAM(jobList, tumorSampleName + "RNA", args.refGenomeFasta, dbSNP = args.dbSNP, intervals = args.intervals, clobber = rnaSplitAndTrim.clobber, emailAddress = args.email, outputDir = outputDir, lastJob = rnaSplitAndTrim, mock = args.mock).returnData
    mPileupRNA = workFlows.MPileup(jobList, tumorSampleName + "RNA", refGenomeFasta = args.refGenomeFasta, gzip = True, emailAddress = args.email, clobber = rnaProcessJobs.clobber, outputDir = outputDir, makeBAMJob = rnaProcessJobs, mock = args.mock).returnData
    #haplotypeCallerRNAJob = workFlows.HaplotypeCaller(jobList, tumorSampleName + "RNA", args.refGenomeFasta, intervals = args.intervals, dbSNP = args.dbSNP, emailAddress = args.email, clobber = rnaProcessJobs.clobber, outputDir = outputDir, lastJob = rnaProcessJobs, mock = args.mock).returnData
    #jointGenotypeRNA = workFlows.JointGenotype(jobList, tumorSampleName + "RNA", args.refGenomeFasta, intervals = args.intervals, dbSNP = args.dbSNP, emailAddress = args.email, standCallConf = 10, clobber = haplotypeCallerRNAJob.clobber, outputDir = outputDir, lastJob = [haplotypeCallerTumorJob, haplotypeCallerRNAJob], mock = args.mock).returnData
    addRNASupportData = workFlows.GetRNASupportMPileup(jobList, tumorSampleName + "RNA", emailAddress = args.email, clobber = mPileupRNA.clobber, outputDir = outputDir, combineSomaticsJob = combineSomaticVariantsJob, rnaMPileupJob = mPileupRNA, mock = args.mock).returnData

    #find any tandem variants (variants at contiguous positions) and combine them into a single polynucleotide variant
    fuseTandemVariants = workFlows.CombineTandemVariants(jobList, tumorSampleName, emailAddress = args.email, clobber = addRNASupportData.clobber, outputDir = outputDir, lastVariantJob = addRNASupportData, mock = args.mock).returnData
    
    #run oncotator and analyze the annotations
    oncotator = workFlows.OncotatorAnalysis(jobList, tumorSampleName, emailAddress = args.email, clobber = fuseTandemVariants.clobber, outputDirectory = outputDir, lastSomaticJob = fuseTandemVariants, mock = args.mock).returnData
    
    #get lists of peptides for each coding variant
    getPeptides = workFlows.GetPeptides(jobList, tumorSampleName, args.peptideDictionaryPickle, sizeList = [8,9,10], emailAddress = args.email, clobber = oncotator.clobber, outputDir = outputDir, oncotatorJob = oncotator, mock = args.mock).returnData

#HLA CALLING SECTION    
    #iteratively set up alignment jobs for normal tissue
    hlaOutputDirectory = outputDir + os.sep + "hlaCalling" + os.sep
    os.mkdir(hlaOutputDirectory)
    hlaAlignJobs = False
    for sample in keys(normalSeqTree): #should be only one here
        for lane in keys(normalSeqTree[sample]):
            currentLane = normalSeqTree[sample][lane]
            if not 2 in currentLane: #which is to say that we did not have a paired end read here
                currentLane[2] = FastqFile(False)
            sampleFile = currentLane[1] #should be a fastqFile object, could use either one, but not guaranteed a second paired end
            if not hlaAlignJobs: #if we have not initiated the job list for this yet, do that
                hlaSickleTrim = workFlows.SickleTrim(jobList, "HLAtrim" + tumorSampleName + str(lane), currentLane[1].filePath, reverseFastqIn = currentLane[2].filePath, emailAddress = args.email, clobber = fuseTandemVariants.clobber, outputDirectory = hlaOutputDirectory, mock = args.mock).returnData
                hlaAlignJobs = workFlows.AlignAlnAddReadGroups(jobList, "HLAaln" + tumorSampleName + str(lane), args.hlaFasta, readGroupSampleName = normalSampleName, readGroupID = lane, clobber = clobber, emailAddress = args.email, lastPe1Job = hlaSickleTrim["forward"], lastPe2Job = hlaSickleTrim["reverse"], flagFilterOut = 4, outputDir = hlaOutputDirectory, mock = args.mock).returnData
                hlaAlignJobs.addJob(workFlows.AlignAlnAddReadGroups(jobList, "HLAalnSing" + tumorSampleName + str(lane), args.hlaFasta, readGroupSampleName = normalSampleName, readGroupID = lane, clobber = clobber, emailAddress = args.email, lastPe1Job = hlaSickleTrim["singleton"], flagFilterOut = 4, outputDir = hlaOutputDirectory, mock = args.mock).returnData)
            else:  #if the job list has been initiated here (so it's not false), add the jobs to it
                hlaSickleTrim = workFlows.SickleTrim(jobList, "HLAtrim" + tumorSampleName + str(lane), currentLane[1].filePath, reverseFastqIn = currentLane[2].filePath, emailAddress = args.email, clobber = fuseTandemVariants.clobber, outputDirectory = hlaOutputDirectory, mock = args.mock).returnData
                hlaAlignJobs.addJob(workFlows.AlignAlnAddReadGroups(jobList, "HLAaln" + tumorSampleName + str(lane), args.hlaFasta, readGroupSampleName = normalSampleName, readGroupID = lane, clobber = clobber, emailAddress = args.email, lastPe1Job = hlaSickleTrim["forward"], lastPe2Job = hlaSickleTrim["reverse"], flagFilterOut = 4, outputDir = hlaOutputDirectory, mock = args.mock).returnData)
                hlaAlignJobs.addJob(workFlows.AlignAlnAddReadGroups(jobList, "HLAalnSing" + tumorSampleName + str(lane), args.hlaFasta, readGroupSampleName = normalSampleName, readGroupID = lane, clobber = clobber, emailAddress = args.email, lastPe1Job = hlaSickleTrim["singleton"], flagFilterOut = 4, outputDir = hlaOutputDirectory, mock = args.mock).returnData)
    #merge and deduplicate
    hlaDeduplicateJobs = workFlows.MergeBAMFilesAndDeduplicate(jobList, "HLA" + tumorSampleName, clobber = hlaAlignJobs.clobber, emailAddress = args.email, outputDir = hlaOutputDirectory, lastJob = hlaAlignJobs, mock = args.mock).returnData
    athlatesJobs = False
    hlaCallTable = {}
    for hlaMolecule in ["A", "B", "C"]:
        if not athlatesJobs:
            currentMoleculeReturn = workFlows.Athlates(jobList, "HLA%s" %hlaMolecule + tumorSampleName,  hlaMolecule, args.athlatesDB, emailAddress = args.email, clobber = hlaDeduplicateJobs.clobber, outputDirectory = outputDir, hlaBAMJob = hlaDeduplicateJobs, mock = args.mock).returnData
            athlatesJobs = currentMoleculeReturn
            hlaCallTable[hlaMolecule] = currentMoleculeReturn.data
        else:
            athlatesJobs = currentMoleculeReturn
            currentMoleculeReturn = workFlows.Athlates(jobList, "HLA%s" %hlaMolecule + tumorSampleName,  hlaMolecule, args.athlatesDB, emailAddress = args.email, clobber = hlaDeduplicateJobs.clobber, outputDirectory = outputDir, hlaBAMJob = hlaDeduplicateJobs, mock = args.mock).returnData
            hlaCallTable[hlaMolecule] = currentMoleculeReturn.data
            athlatesJobs.addJob(currentMoleculeReturn)
            
#Integrate the HLA data and the peptide data to submit to netMHC
    addHLA = workFlows.AddHLA(jobList, tumorSampleName, hlaCallTable, emailAddress = args.email, clobber = currentMoleculeReturn.clobber, outputDir = outputDir, lastSomaticPickleJob = getPeptides, hlaJobs = athlatesJobs, mock = args.mock).returnData
    
#Run netMHC predictions on the peptides and create an output file
    netMHC = workFlows.NetMHCAnalysis(jobList, tumorSampleName, emailAddress = args.email, clobber = addHLA.clobber, outputDir = outputDir, lastSomaticPickleJob = addHLA, mock = args.mock).returnData
  
#CLOSING OUT    
    #Capstone to mark completion
    capstoneDependencies = [netMHC.jobID]
    capstoneJobID = workFlows.Capstone(jobList, tumorSampleName, capstoneDependencies, args.email, outputDir, args.mock)
    
    #ending the pipeline
    pipeline.end()

if __name__ == '__main__':
    main()

