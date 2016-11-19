#!/usr/bin/env python3

import os
global gatkPath
gatkPath = os.getcwd() + "/bin/GenomeAnalysisTK.jar"

class IndelRealignment(object):
    
    def __init__(self, sampleName, bamIn, refGenomeFasta, known = False, bedFile = False, allowPotentiallyMisencodedQualityScores = True, clobber = False, outputDirectory = ""):
        import runnerSupport
        self.sampleName = sampleName
        self.allowPotentiallyMisencodedQualityScores = allowPotentiallyMisencodedQualityScores
        self.refGenomeFasta = refGenomeFasta
        if type(known) == str:
            known = [known]
        self.known = known
        self.bedFile = bedFile
        if not outputDirectory:
            self.outputDirectory = ""
        else:
            if not os.path.isdir(outputDirectory):
                raise RuntimeError("Output directory %s does not exist.  Please make the directory before creating jobs." %(outputDirectory))
            if not outputDirectory.endswith(os.sep):
                self.outputDirectory = outputDirectory + os.sep
            else:
                self.outputDirectory = outputDirectory
        self.bamIn = bamIn
        #SANITY TEST ALL THE THINGS
        runnerSupport.checkForRequiredFile(self.bamIn, "BAM file for realignment")
        if self.known:
            for file in self.known:
                runnerSupport.checkForRequiredFile(file, "Known variant reference VCF " + file)
        if self.bedFile:
            runnerSupport.checkForRequiredFile(self.bedFile, "BED file containing intervals")
        self.checkRefGenome()
        #DONE SANITY CHECKING. FOR NOW.
        self.makeAndCheckOutputFileNames()
        self.targetCreationCommand, self.realignerCommand = self.createGATKCommand()
        
    def makeAndCheckOutputFileNames(self):
        import runnerSupport
        self.intervalsOut = self.outputDirectory + runnerSupport.stripDirectoryAndExtension(self.bamIn) + ".intervals"
        self.bamOut = self.outputDirectory + runnerSupport.stripDirectoryAndExtension(self.bamIn) + ".indelRealign" + ".bam"
        self.clobber = runnerSupport.checkForOverwriteRisk(self.bamOut, self.sampleName, self.clobber)
        self.clobber = runnerSupport.checkForOverwriteRisk(self.intervalsOut, self.sampleName, self.clobber)
        
    def checkRefGenome(self):
        import os
        import runnerSupport
        runnerSupport.checkForRequiredFile(self.refGenomeFasta, "reference genome file")
        runnerSupport.checkForRequiredFile(self.refGenomeFasta + ".fai", "reference genome fasta index", "Please move the index to this location or use samtools faidx to create one.")
        runnerSupport.checkForRequiredFile(".".join(self.refGenomeFasta.split(".")[:-1]) + ".dict", "reference genome dictionary", "Please move the index to this location or use picard create sequence dictionary to create one.")

    def createGATKCommand(self):
        import runnerSupport
        flagValuesTarget = {"-T" : "RealignerTargetCreator",
                            "-R" : self.refGenomeFasta,
                            "-I" : self.bamIn,
                            "-o" : self.intervalsOut,
                            "-L" : self.bedFile,
                            "flaggedlist" : ["-known", self.known],
                            "--allow_potentially_misencoded_quality_scores" : self.allowPotentiallyMisencodedQualityScores}
        flagValuesRealign = {"-T" : "IndelRealigner",
                             "-R" : self.refGenomeFasta,
                             "-I" : self.bamIn,
                             "-o" : self.bamOut,
                             "-L" : self.bedFile,
                             "-targetIntervals" : self.intervalsOut,
                             "flaggedlist" : ["-known", self.known],
                             "--allow_potentially_misencoded_quality_scores" : self.allowPotentiallyMisencodedQualityScores}
        gatkArgs = [programPaths["java"], "-Xmx1g", "-jar", gatkPath, flagValuesTarget]
        argumentFormatter = runnerSupport.ArgumentFormatter(gatkArgs)
        targetCommand = argumentFormatter.argumentString
        gatkArgs = [programPaths["java"], "-Xmx1g", "-jar", gatkPath, flagValuesRealign]
        argumentFormatter = runnerSupport.ArgumentFormatter(gatkArgs)
        realignCommand = argumentFormatter.argumentString
        return (targetCommand, realignCommand)
        
class BQSR(object):
    
    def __init__(self, sampleName, bamIn, refGenomeFasta, known, bedFile = False, allowPotentiallyMisencodedQualityScores = True, clobber = False, outputDirectory = ""):
        import runnerSupport
        self.sampleName = sampleName
        self.allowPotentiallyMisencodedQualityScores = allowPotentiallyMisencodedQualityScores
        self.refGenomeFasta = refGenomeFasta
        if type(known) == str:
            known = [known]
        self.known = known
        if not outputDirectory:
            self.outputDirectory = ""
        else:
            if not os.path.isdir(outputDirectory):
                raise RuntimeError("Output directory %s does not exist.  Please make the directory before creating jobs." %(outputDirectory))
            if not outputDirectory.endswith(os.sep):
                self.outputDirectory = outputDirectory + os.sep
            else:
                self.outputDirectory = outputDirectory
        self.bamIn = bamIn
        #SANITY TEST ALL THE THINGS
        runnerSupport.checkForRequiredFile(self.bamIn, "BAM file for realignment")
        if self.known:
            for file in self.known:
                runnerSupport.checkForRequiredFile(file, "Known variant reference VCF " + file)
        if self.bedFile:
            runnerSupport.checkForRequiredFile(self.bedFile, "BED file containing target intervals")
        self.checkRefGenome()
        #DONE SANITY CHECKING. FOR NOW.
        self.makeAndCheckOutputFileNames()
        self.firstPassCommand, self.secondPassCommand, self.plotCommand, self.printReadsCommand = self.createGATKCommand()
        
    def makeAndCheckOutputFileNames(self):
        import runnerSupport
        self.recalTableOut = self.outputDirectory + runnerSupport.stripDirectoryAndExtension(self.bamIn) + ".recal.table"
        self.postRecalTableOut = self.outputDirectory + runnerSupport.stripDirectoryAndExtension(self.bamIn) + ".post_recal.table"
        self.plotsOut = self.outputDirectory + runnerSupport.stripDirectoryAndExtension(self.bamIn) + ".base_quality_recalibration_plots.pdf"
        self.bamOut = self.outputDirectory + runnerSupport.stripDirectoryAndExtension(self.bamIn) + ".recal" + ".bam"
        self.clobber = runnerSupport.checkForOverwriteRisk(self.bamOut, self.sampleName, self.clobber)
        self.clobber = runnerSupport.checkForOverwriteRisk(self.recalTableOut, self.sampleName, self.clobber)
        self.clobber = runnerSupport.checkForOverwriteRisk(self.postRecalTableOut, self.sampleName, self.clobber)
        self.clobber = runnerSupport.checkForOverwriteRisk(self.plotsOut, self.sampleName, self.clobber)
        
    def checkRefGenome(self):
        import os
        runnerSupport.checkForRequiredFile(self.refGenomeFasta, "reference genome file")
        runnerSupport.checkForRequiredFile(self.refGenomeFasta + ".fai", "reference genome fasta index", "Please move the index to this location or use samtools faidx to create one.")
        runnerSupport.checkForRequiredFile(".".join(self.refGenomeFasta.split(".")[:-1]) + ".dict", "reference genome dictionary", "Please move the index to this location or use picard create sequence dictionary to create one.")

    def createGATKCommand(self):
        import runnerSupport
        flagValuesPre =    {"-T" : "BaseRecalibrator",
                            "-R" : self.refGenomeFasta,
                            "-I" : self.bamIn,
                            "-o" : self.recalTableOut,
                            "-L" : self.bedFile,
                            "flaggedlist" : ["-known", self.known],
                            "--allow_potentially_misencoded_quality_scores" : self.allowPotentiallyMisencodedQualityScores}
        flagValuesPost =   {"-T" : "BaseRecalibrator",
                            "-R" : self.refGenomeFasta,
                            "-I" : self.bamIn,
                            "-o" : self.postRecalTableOut,
                            "-L" : self.bedFile,
                            "-BQSR" : self.recalTableOut,
                            "flaggedlist" : ["-known", self.known],
                            "--allow_potentially_misencoded_quality_scores" : self.allowPotentiallyMisencodedQualityScores}
        flagValuesPlot =   {"-T" : "AnalyzeCovariates",
                            "-R" : self.refGenomeFasta,
                            "-I" : self.bamIn,
                            "-plots" : self.plotsOut,
                            "-L" : self.bedFile,
                            "-before" : self.recalTableOut,
                            "-after" : self.postRecalTableOut,
                            "--allow_potentially_misencoded_quality_scores" : self.allowPotentiallyMisencodedQualityScores}
        flagValuesPrint =  {"-T" : "PrintReads",
                            "-R" : self.refGenomeFasta,
                            "-I" : self.bamIn,
                            "-o" : self.bamOut,
                            "-L" : self.bedFile,
                            "-BQSR" : self.recalTableOut,
                            "--allow_potentially_misencoded_quality_scores" : self.allowPotentiallyMisencodedQualityScores}
        gatkArgs = [programPaths["java"], "-Xmx1g", "-jar", gatkPath, flagValuesPre]
        argumentFormatter = runnerSupport.ArgumentFormatter(gatkArgs)
        firstPass = argumentFormatter.argumentString
        gatkArgs = [programPaths["java"], "-Xmx1g", "-jar", gatkPath, flagValuesPost]
        argumentFormatter = runnerSupport.ArgumentFormatter(gatkArgs)
        secondPass = argumentFormatter.argumentString
        gatkArgs = [programPaths["java"], "-Xmx1g", "-jar", gatkPath, flagValuesPlot]
        argumentFormatter = runnerSupport.ArgumentFormatter(gatkArgs)
        plot = argumentFormatter.argumentString
        gatkArgs = [programPaths["java"], "-Xmx1g", "-jar", gatkPath, flagValuesPrint]
        argumentFormatter = runnerSupport.ArgumentFormatter(gatkArgs)
        printReads = argumentFormatter.argumentString
        return (firstPass, secondPass, plot, printReads)
    
class HaplotypeCaller(object):
    
    def __init__(self, sampleName, bamIn, refGenomeFasta, bedFile = False, emitRefConfidence = "GVCF", variantIndexType = "Linear", variantIndexParameter = 128000, dbSNP = False, annotation = ["Coverage", "AlleleBalance"], intervalPadding = 100, pairHiddenMarkovModel = "VECTOR_LOGLESS_CACHING", allowPotentiallyMisencodedQualityScores = True, clobber = False, outputDirectory = ""):
        import runnerSupport
        self.bedFile = bedFile
        self.sampleName = sampleName
        self.allowPotentiallyMisencodedQualityScores = allowPotentiallyMisencodedQualityScores
        self.refGenomeFasta = refGenomeFasta
        self.emitRefConfidence = emitRefConfidence
        self.variantIndexType = variantIndexType
        self.variantIndexParameter = variantIndexParameter
        self.dbSNP = dbSNP
        self.annotation = annotation
        if not type(self.annotation) == list:
            self.annotation = [self.annotation]
        self.intervalPadding = intervalPadding
        self.pairHiddenMarkovModel = pairHiddenMarkovModel
        if not outputDirectory:
            self.outputDirectory = ""
        else:
            if not os.path.isdir(outputDirectory):
                raise RuntimeError("Output directory %s does not exist.  Please make the directory before creating jobs." %(outputDirectory))
            if not outputDirectory.endswith(os.sep):
                self.outputDirectory = outputDirectory + os.sep
            else:
                self.outputDirectory = outputDirectory
        self.bamIn = bamIn
        #SANITY TEST ALL THE THINGS
        runnerSupport.checkForRequiredFile(self.bamIn, "BAM file for realignment")
        if self.dbSNP:
            runnerSupport.checkForRequiredFile(dbSNP, "dbSNP file")
        if self.bedFile:
            runnerSupport.checkForRequiredFile(self.bedFile, "BED file containing target intervals")
        if not type(self.emitRefConfidence) == int and not self.emitRefConfidence == "GVCF":
            raise RuntimeError("Invalid emitRefConfidence value passed to haplotype caller: %s" %(self.emitRefConfidence))
        if not type(self.variantIndexParameter):
            raise RuntimeError("Variant Index Parameter value was not an integer.")
        if not type(self.intervalPadding) == int:
            raise RuntimeError("Interval padding value for haplotype caller must be an integer.")
        self.checkRefGenome()
        #DONE SANITY CHECKING. FOR NOW.
        self.makeAndCheckOutputFileNames()
        self.haplotypeCallerCommand = self.createGATKCommand()
        
    def makeAndCheckOutputFileNames(self):
        import runnerSupport
        self.gvcfOut = self.outputDirectory + self.sampleName + ".gvcf"
        self.clobber = runnerSupport.checkForOverwriteRisk(self.gvcfOut, self.sampleName, self.clobber)
        
    def checkRefGenome(self):
        import os
        import runnerSupport
        runnerSupport.checkForRequiredFile(self.refGenomeFasta, "reference genome file")
        runnerSupport.checkForRequiredFile(self.refGenomeFasta + ".fai", "reference genome fasta index", "Please move the index to this location or use samtools faidx to create one.")
        runnerSupport.checkForRequiredFile(".".join(self.refGenomeFasta.split(".")[:-1]) + ".dict", "reference genome dictionary", "Please move the index to this location or use picard create sequence dictionary to create one.")

    def createGATKCommand(self):
        import runnerSupport
        flagValues = {"-T" : "HaplotypeCaller",
                      "-R" : self.refGenomeFasta,
                      "-I" : self.bamIn,
                      "-o" : self.gvcfOut,
                      "-L" : self.bedFile,
                      "--variant_index_type" : self.variantIndexType,
                      "--variant_index_parameter" : self.variantIndexParameter,
                      "-D" : self.dbSNP,
                      "--interval_padding" : self.intervalPadding,
                      "-pairHMM" : self.pairHiddenMarkovModel,
                      "flaggedlist" : ["--annotation", self.annotation],
                      "--allow_potentially_misencoded_quality_scores" : self.allowPotentiallyMisencodedQualityScores}
        gatkArgs = [programPaths["java"], "-Xmx1g", "-jar", gatkPath, flagValues]
        argumentFormatter = runnerSupport.ArgumentFormatter(gatkArgs)
        haplotypeCallerCommand = argumentFormatter.argumentString
        return (haplotypeCallerCommand)
    
class DepthOfCoverage(object):
    
    def __init__(self, sampleName, bamIn, refGenomeFasta, bedFile = False, allowPotentiallyMisencodedQualityScores = True, clobber = False, outputDirectory = ""):
        import runnerSupport
        self.bedFile = bedFile
        self.sampleName = sampleName
        self.allowPotentiallyMisencodedQualityScores = allowPotentiallyMisencodedQualityScores
        self.refGenomeFasta = refGenomeFasta
        if not outputDirectory:
            self.outputDirectory = ""
        else:
            if not os.path.isdir(outputDirectory):
                raise RuntimeError("Output directory %s does not exist.  Please make the directory before creating jobs." %(outputDirectory))
            if not outputDirectory.endswith(os.sep):
                self.outputDirectory = outputDirectory + os.sep
            else:
                self.outputDirectory = outputDirectory
        self.bamIn = bamIn
        #SANITY TEST ALL THE THINGS
        runnerSupport.checkForRequiredFile(self.bamIn, "BAM file for depth check")
        if self.bedFile:
            runnerSupport.checkForRequiredFile(self.bedFile, "BED file containing target intervals")
        self.checkRefGenome()
        #DONE SANITY CHECKING. FOR NOW.
        self.makeAndCheckOutputFileNames()
        self.depthOfCoverageCommand = self.createGATKCommand()
        
    def makeAndCheckOutputFileNames(self):
        self.depthOut = self.outputDirectory + self.sampleName + ".depthCov"
        self.clobber = runnerSupport.checkForOverwriteRisk(self.depthOut, self.sampleName, self.clobber)
        
    def checkRefGenome(self):
        import os
        import runnerSupport
        runnerSupport.checkForRequiredFile(self.refGenomeFasta, "reference genome file")
        runnerSupport.checkForRequiredFile(self.refGenomeFasta + ".fai", "reference genome fasta index", "Please move the index to this location or use samtools faidx to create one.")
        runnerSupport.checkForRequiredFile(".".join(self.refGenomeFasta.split(".")[:-1]) + ".dict", "reference genome dictionary", "Please move the index to this location or use picard create sequence dictionary to create one.")

    def createGATKCommand(self):
        import runnerSupport
        flagValues = {"-T" : "DepthOfCoverage",
                      "-R" : self.refGenomeFasta,
                      "-I" : self.bamIn,
                      "-o" : self.depthOut,
                      "-L" : self.bedFile,
                      "--allow_potentially_misencoded_quality_scores" : self.allowPotentiallyMisencodedQualityScores}
        gatkArgs = [programPaths["java"], "-Xmx1g", "-jar", gatkPath, flagValues]
        argumentFormatter = runnerSupport.ArgumentFormatter(gatkArgs)
        depthCommand = argumentFormatter.argumentString
        return (depthCommand)
    
class Mutect1(object):
    
    def __init__(self, sampleName, normalBamIn, tumorBamIn, refGenomeFasta, bedFile = False, dbSNP = False, allowPotentiallyMisencodedQualityScores = True, clobber = False, outputDirectory = ""):
        import runnerSupport
        self.sampleName = sampleName
        self.allowPotentiallyMisencodedQualityScores = allowPotentiallyMisencodedQualityScores
        self.refGenomeFasta = refGenomeFasta
        self.normalBamIn = normalBamIn
        self.tumorBamIn = tumorBamIn
        self.bedFile = bedFile
        self.dbSNP = dbSNP
        if not outputDirectory:
            self.outputDirectory = ""
        else:
            if not os.path.isdir(outputDirectory):
                raise RuntimeError("Output directory %s does not exist.  Please make the directory before creating jobs." %(outputDirectory))
            if not outputDirectory.endswith(os.sep):
                self.outputDirectory = outputDirectory + os.sep
            else:
                self.outputDirectory = outputDirectory
        #SANITY TEST ALL THE THINGS
        runnerSupport.checkForRequiredFile(self.normalBamIn, "Normal BAM for comparison")
        runnerSupport.checkForRequiredFile(self.tumorBamIn, "Tumor BAM for comparison")
        if self.dbSNP:
            runnerSupport.checkForRequiredFile(dbSNP, "dbSNP file")
        if self.bedFile:
            runnerSupport.checkForRequiredFile(self.bedFile, "BED file containing target intervals")
        runnerSupport.checkTypes([self.allowPotentiallyMisencodedQualityScores], bool)
        self.checkRefGenome()
        #DONE SANITY CHECKING. FOR NOW.
        self.makeAndCheckOutputFileNames()
        self.mutect1Command = self.createGATKCommand()
        
    def makeAndCheckOutputFileNames(self):
        import runnerSupport
        self.vcfOut = self.outputDirectory + self.sampleName + ".mutect.vcf"
        self.clobber = runnerSupport.checkForOverwriteRisk(self.vcfOut, self.sampleName, self.clobber)
        
    def checkRefGenome(self):
        import os
        import runnerSupport
        runnerSupport.checkForRequiredFile(self.refGenomeFasta, "reference genome file")
        runnerSupport.checkForRequiredFile(self.refGenomeFasta + ".fai", "reference genome fasta index", "Please move the index to this location or use samtools faidx to create one.")
        runnerSupport.checkForRequiredFile(".".join(self.refGenomeFasta.split(".")[:-1]) + ".dict", "reference genome dictionary", "Please move the index to this location or use picard create sequence dictionary to create one.")

    def createGATKCommand(self):
        import runnerSupport
        flagValues = {"-T" : "Mutect",
                      "-R" : self.refGenomeFasta,
                      "--input_file:normal" : self.normalBamIn,
                      "--input_file:tumor" : self.tumorBamIn,
                      "-o" : self.vcfOut,
                      "-L" : self.bedFile,
                      "-D" : self.dbSNP,
                      "--allow_potentially_misencoded_quality_scores" : self.allowPotentiallyMisencodedQualityScores}
        gatkArgs = [programPaths["java"], "-Xmx1g", "-jar", gatkPath, flagValues]
        argumentFormatter = runnerSupport.ArgumentFormatter(gatkArgs)
        mutect1Command = argumentFormatter.argumentString
        return (mutect1Command)
    
class JointGenotype(object):
    
    def __init__(self, sampleName, gvcfIn, comparisonGVCFDir, refGenomeFasta, dbSNP = False,  allowPotentiallyMisencodedQualityScores = True, clobber = False, outputDirectory = ""):
        import runnerSupport
        self.sampleName = sampleName
        self.allowPotentiallyMisencodedQualityScores = allowPotentiallyMisencodedQualityScores
        self.refGenomeFasta = refGenomeFasta
        self.comparisonGVCFDir = comparisonGVCFDir
        self.gvcfIn = gvcfIn
        self.dbSNP = dbSNP
        if not outputDirectory:
            self.outputDirectory = ""
        else:
            if not os.path.isdir(outputDirectory):
                raise RuntimeError("Output directory %s does not exist.  Please make the directory before creating jobs." %(outputDirectory))
            if not outputDirectory.endswith(os.sep):
                self.outputDirectory = outputDirectory + os.sep
            else:
                self.outputDirectory = outputDirectory
        #SANITY TEST ALL THE THINGS
        runnerSupport.checkForRequiredFile(self.gvcfIn, "GVCF for joint genotype calling")
        if self.dbSNP:
            runnerSupport.checkForRequiredFile(dbSNP, "dbSNP file")
        runnerSupport.checkTypes([self.allowPotentiallyMisencodedQualityScores], bool)
        self.checkRefGenome()
        #DONE SANITY CHECKING. FOR NOW.
        self.makeGVCFList()
        self.makeAndCheckOutputFileNames()
        self.jointGenotypeCommand = self.createGATKCommand()
        
    def makeGVCFList(self):
        import os
        if not self.comparisonGVCFDir.endswith(os.sep):
            self.comparisonGVCFDir += os.sep
        if not os.path.isdir(self.comparisonGVCFDir):
            raise RuntimeError("GVCF comparison set directory not found at %s" %(self.comparisonGVCFDir))
        rawDirectoryListing = os.listdir(self.comparisonGVCFDir)
        self.gvcfList = []
        for file in rawDirectoryListing:
            if os.path.isfile(self.comparisonGVCFDir + file):
                if file.lower().endswith("g.vcf") or file.lower().endswith(".gvcf"):
                    self.gvcfList.append(self.comparisonGVCFDir + file)
        if not self.gvcfList:
            raise RuntimeError("No GVCFs were present in the specified comparison GVCF directory. Joint genotyping is impossible. Directory: %s" %(self.comparisonGVCFDir))
                
    def makeAndCheckOutputFileNames(self):
        import runnerSupport
        self.vcfOut = self.outputDirectory + self.sampleName + ".vcf"
        self.clobber = runnerSupport.checkForOverwriteRisk(self.vcfOut, self.sampleName, self.clobber)
        
    def checkRefGenome(self):
        import os
        import runnerSupport
        runnerSupport.checkForRequiredFile(self.refGenomeFasta, "reference genome file")
        runnerSupport.checkForRequiredFile(self.refGenomeFasta + ".fai", "reference genome fasta index", "Please move the index to this location or use samtools faidx to create one.")
        runnerSupport.checkForRequiredFile(".".join(self.refGenomeFasta.split(".")[:-1]) + ".dict", "reference genome dictionary", "Please move the index to this location or use picard create sequence dictionary to create one.")

    def createGATKCommand(self):
        import runnerSupport
        flagValues = {"-T" : "GenotypeGVCFs",
                      "-R" : self.refGenomeFasta,
                      "--variant" : self.gvcfIn,
                      "flaggedlist" : ["--variant", self.gvcfList],
                      "-o" : self.vcfOut,
                      "-D" : self.dbSNP,
                      "--allow_potentially_misencoded_quality_scores" : self.allowPotentiallyMisencodedQualityScores}
        gatkArgs = [programPaths["java"], "-Xmx1g", "-jar", gatkPath, flagValues]
        argumentFormatter = runnerSupport.ArgumentFormatter(gatkArgs)
        jointGenotypeCommand = argumentFormatter.argumentString
        return (jointGenotypeCommand)

class VQSRResource(object):
    
    def __init__(self, fileName, resourceName, known = False, training = False, truth = False, prior = False):
        import runnerSupport
        self.fileName = fileName
        self.resourceName = resourceName
        self.known = known
        self.training = training
        self.truth = truth
        self.prior = prior
        runnerSupport.checkTypes([known, training, truth], bool)
        runnerSupport.checkTypes(prior, float)
        runnerSupport.checkTypes(resourceName, str)
        runnerSupport.checkForRequiredFile(fileName, "Reference file for VQSR resource %s" %(self.resourceName))
        self.argument = self.formArgument()
        
    def boolPrint(self, value):
        if value:
            return "true"
        else:
            return "false"
        
    def formArgument(self):
        known = "=".join("known", self.boolPrint(self.known))
        training = "=".join("training", self.boolPrint(self.training))
        truth = "=".join("truth", self.boolPrint(self.truth))
        prior = "=".join("prior", str(self.prior))
        descriptors = ",".join(self.resourceName, known, training, truth)
        return descriptors + " " + self.fileName
    
    def __str__(self):
        return self.argument
    
class VQSR(object):
    
    def __init__(self, sampleName, vcfIn, resourceList, refGenomeFasta, maxGaussians = "default", annotationsList = False, mode = False, truthSensitivityFilter = False, allowPotentiallyMisencodedQualityScores = True, clobber = False, outputDirectory = ""):
        import runnerSupport
        snpAnnotationStandard = ["QD", "MQ", "MQRankSum", "ReadPosRankSum", "FS", "SOR"]
        indelAnnotationStandard = ["QD", "MQRankSum", "ReadPosRankSum", "FS", "SOR"]
        snpDefaultMaxGaussians = False
        indelDefaultMaxGaussians = 4
        snpDefaultTruthSensitivityFilter = 99.5
        indelDefaultTruthSensitivityFilter = 99.0
        self.mode = mode
        if not mode:
            raise RuntimeError("VQSR mode must be specified.  None was.")
        if maxGaussians == "default":
            if self.mode == "SNP":
                self.maxGaussians = snpDefaultMaxGaussians
            elif self.mode == "INDEL":
                self.maxGaussians = indelDefaultMaxGaussians
        else:
            self.maxGaussians = maxGaussians
        if not truthSensitivityFilter:
            if self.mode == "SNP":
                self.truthSensitivityFilter = snpDefaultTruthSensitivityFilter
            elif self.mode == "INDEL":
                self.truthSensitivityFilter = indelDefaultTruthSensitivityFilter
        else:
            self.truthSensitivityFilter = truthSensitivityFilter
        runnerSupport.checkTypes(mode, str)
        self.mode = self.mode.upper()
        if not self.mode in ["INDEL", "SNP"]:
            raise RuntimeError("VQSR mode must be set to SNP or INDEL only.")   #validating mode so we can check on annotations list
        if annotationsList:
            self.annotationsList = annotationsList
        else:
            if self.mode == "SNP":
                self.annotationsList = snpAnnotationStandard
            elif self.mode == "INDEL":
                self.annotationsList = indelAnnotationStandard
        self.sampleName = sampleName
        self.vcfIn = vcfIn
        self.resourceList = resourceList
        self.allowPotentiallyMisencodedQualityScores = allowPotentiallyMisencodedQualityScores
        self.refGenomeFasta = refGenomeFasta
        if not outputDirectory:
            self.outputDirectory = ""
        else:
            if not os.path.isdir(outputDirectory):
                raise RuntimeError("Output directory %s does not exist.  Please make the directory before creating jobs." %(outputDirectory))
            if not outputDirectory.endswith(os.sep):
                self.outputDirectory = outputDirectory + os.sep
            else:
                self.outputDirectory = outputDirectory
        #SANITY TEST ALL THE THINGS
        runnerSupport.checkForRequiredFile(self.vcfIn, "VCF for analysis")
        runnerSupport.checkTypes(annotationsList, list)
        runnerSupport.checkTypes(truthSensitivityFilter, float)
        for resource in self.resourceList:
            runnerSupport.checkTypes(resource, VQSRResource)
        self.checkRefGenome()
        #DONE SANITY CHECKING. FOR NOW.
        self.makeAndCheckOutputFileNames()
        self.analysisCommand, self.executionCommand = self.createGATKCommand()
        
    def makeAndCheckOutputFileNames(self):
        import runnerSupport
        #checking only the vcf output and not any of the secondary outputs.  If we are creating a new VCF, we need to create new secondary outputs to go with it and SHOULD overwrite to avoid confusion.
        newExtension = "." + self.mode.lower() + "Recal.vcf"
        self.vcfOut = self.outputDirectory + runnerSupport.stripDirectoryAndExtension(self.vcfIn) + newExtension
        self.recalFile = ".".join(self.sampleName, self.mode, "recal")
        self.tranchesFile = ".".join(self.sampleName, self.mode, "tranches")
        self.rScriptFile = ".".join(self.sampleName, self.mode, "plots.R")
        self.clobber = runnerSupport.checkForOverwriteRisk(self.vcfOut, self.sampleName, self.clobber)
        
    def checkRefGenome(self):
        import os
        import runnerSupport
        runnerSupport.checkForRequiredFile(self.refGenomeFasta, "reference genome file")
        runnerSupport.checkForRequiredFile(self.refGenomeFasta + ".fai", "reference genome fasta index", "Please move the index to this location or use samtools faidx to create one.")
        runnerSupport.checkForRequiredFile(".".join(self.refGenomeFasta.split(".")[:-1]) + ".dict", "reference genome dictionary", "Please move the index to this location or use picard create sequence dictionary to create one.")

    def createGATKCommand(self):
        import runnerSupport
        resourceArguments = [resource.argument for resource in self.resourceList]
        flagValuesAnalyze = {"-T" : "VariantRecalibrator",
                             "-R" : self.refGenomeFasta,
                             "-input" : self.vcfIn,
                             "flaggedlist1" : ["-resource", resourceArguments, ":"],
                             "flaggedlist2" : ["-an", self.annotationsList],
                             "-mode" : self.mode,
                             "-maxGaussians" : self.maxGaussians,
                             "-recalFile" : self.recalFile,
                             "-tranchesFile" : self.tranchesFile,
                             "rscriptFile" : self.rScriptFile,
                             "--allow_potentially_misencoded_quality_scores" : self.allowPotentiallyMisencodedQualityScores}
        flagValuesApply  =  {"-T" : "ApplyRecalibration",
                             "-R" : self.refGenomeFasta,
                             "-input" : self.vcfIn,
                             "-mode" : self.mode,
                             "-recalFile" : self.recalFile,
                             "-tranchesFile" : self.tranchesFile,
                             "--out" : self.vcfOut,
                             "--ts_filter_level" : self.truthSensitivityFilter,
                             "--allow_potentially_misencoded_quality_scores" : self.allowPotentiallyMisencodedQualityScores}
        gatkArgs = [programPaths["java"], "-Xmx1g", "-jar", gatkPath, flagValuesAnalyze]
        argumentFormatter = runnerSupport.ArgumentFormatter(gatkArgs)
        analysisCommand = argumentFormatter.argumentString
        gatkArgs = [programPaths["java"], "-Xmx1g", "-jar", gatkPath, flagValuesApply]
        argumentFormatter = runnerSupport.ArgumentFormatter(gatkArgs)
        executionCommand = argumentFormatter.argumentString
        return (analysisCommand, executionCommand)
