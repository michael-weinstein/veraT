#!/usr/bin/env python3

import os
global programPaths
programPaths = {"bwa" : os.getcwd() + "/bin/bwa-0.7.15/bwa",
                "java" : os.getcwd() + "/bin/jre1.8.0_77/bin/java",
                "samtools" : "/u/local/apps/samtools/1.2/gcc-4.4.7/bin/samtools",
                "extractVariants" : os.getcwd() + "/analysisScripts/extractVariants.py",
                "combineVariants" : os.getcwd() + "/analysisScripts/combineExtractedVariants.py",
                "python3" : "/u/local/apps/python/3.4.3/bin/python3",                
                "bgzip" : os.getcwd() + "/bin/tabix/tabix-0.2.6/bgzip",
                "tabix" : os.getcwd() + "/bin/tabix/tabix-0.2.6/tabix",
                "varscan" : os.getcwd() + "/bin/VarScan.v2.4.0.jar",
                "bam-readcount" : os.getcwd() + "/bin/bamReadCount/bin/bam-readcount",
                "vcfReader" : os.getcwd() + "/runners/variantReaders/vcfReader.py",
                "mutectReader" : os.getcwd() + "/runners/variantReaders/mutectReader.py",
                "varScanPositions" : os.getcwd() + "/runners/variantReaders/varScanPositionPuller.py",
                "varScanReader" : os.getcwd() + "/runners/variantReaders/varScanReader.py",
                "perl" : "/usr/bin/perl",
                "varScanFPFilter" : os.getcwd() + "/runners/variantReaders/fpfilter.pl",
                "variantCombine" : os.getcwd() + "/runners/variantReaders/variantCombine.py",
                "getRNASupport" : os.getcwd() + "/runners/variantReaders/getRNASupport.py",
                "hisat2" : os.getcwd() + "/bin/hisat2/hisat2"}

class BWAlign(object):
    
    def __init__(self, sampleName, refGenomeFasta, pairedEnd1File, pairedEnd2File = "", cores = 1, qualityThreshold = 5, seedLength = 32, maxMismatchesPerSeed = 2, maximumGapOpens = 1, clobber = None, outputDirectory = ""):
        import os
        import runnerSupport
        if not outputDirectory:
            self.outputDirectory = ""
        else:
            if not os.path.isdir(outputDirectory):
                raise RuntimeError("Output directory %s does not exist.  Please make the directory before creating jobs." %(outputDirectory))
            if not outputDirectory.endswith(os.sep):
                self.outputDirectory = outputDirectory + os.sep
            else:
                self.outputDirectory = outputDirectory
        self.sampleName = sampleName
        self.pe1 = pairedEnd1File
        self.pe2 = pairedEnd2File
        self.cores = cores
        self.threads = cores
        self.qualityThreshold = qualityThreshold
        self.seedLength = seedLength
        self.maxMismatchesPerSeed = maxMismatchesPerSeed
        self.maximumGapOpens = maximumGapOpens
        self.refGenomeFasta = refGenomeFasta
        self.clobber = clobber
        #SANITY CHECK ALL THE THINGS
        if not type(self.pe1) == str or not self.pe1:
            raise RuntimeError("Paired end 1 file must be specified and must be a string.  Value passed was \"%s\"" %(self.pe1))
        runnerSupport.checkForRequiredFile(self.pe1, "paired end 1 file")
        if not type(self.pe2) == str:
            raise RuntimeError("Paired end 2 file must be specified and must be a string.  Value passed was \"%s\"" %(self.pe2))
        if self.pe2:
            runnerSupport.checkForRequiredFile(self.pe2, "paired end 2 file")
            self.pairedEnd = True
        else:
            self.pairedEnd = False
        if not type(self.threads) == int and self.threads > 0:
            raise RuntimeError("Thread count must be a positive integer. %s was passed." %(self.threads))
        if not type(self.qualityThreshold) == int and self.qualityThreshold > 0:
            raise RuntimeError("Thread count must be a positive integer. %s was passed." %(self.qualityThreshold))
        if not type(self.seedLength) == int and self.seedLength > 0:
            raise RuntimeError("Thread count must be a positive integer. %s was passed." %(self.seedLength))
        if not type(self.maxMismatchesPerSeed) == int and self.maxMismatchesPerSeed > 0:
            raise RuntimeError("Thread count must be a positive integer. %s was passed." %(self.maxMismatchesPerSeed))
        if not type(self.maximumGapOpens) == int and self.maximumGapOpens > 0:
            raise RuntimeError("Thread count must be a positive integer. %s was passed." %(self.maximumGapOpens))
        if not type(self.refGenomeFasta) == str or not self.refGenomeFasta:
            raise RuntimeError("Reference genome fasta must be passed as a string. %s was passed." %(self.refGenomeFasta))
        self.checkRefGenome()
        self.makeAndCheckOutputFileNames()
        #DONE SANITY CHECKING ALL THE THINGS. FOR NOW.
        self.pe1Command = self.makeBWAlignCommand(self.pe1, self.pe1Out)
        if not self.pairedEnd:
            self.pe2Command = ""
        else:
            self.pe2Command = self.makeBWAlignCommand(self.pe2, self.pe2Out)
        self.samCommand = self.makeSAMCommand()
    
    def checkRefGenome(self):
        import os
        import runnerSupport
        runnerSupport.checkForRequiredFile(self.refGenomeFasta, "reference genome file")
        runnerSupport.checkForRequiredFile(self.refGenomeFasta + ".fai", "reference genome fasta index", "Please move the index to this location or use samtools faidx to create one.")
        bwaIndexFileExtensions = [".amb", ".ann", ".bwt", ".pac", ".sa"]
        for extension in bwaIndexFileExtensions:
            indexFile = self.refGenomeFasta + extension
            runnerSupport.checkForRequiredFile(indexFile, "one or more BWA index files", "Please move the index to this location or use bwa index to (re)create one.")
    
    def makeAndCheckOutputFileNames(self):
        import runnerSupport
        self.pe1Out = self.outputDirectory + self.sampleName + ".pe1.sai"
        self.clobber = runnerSupport.checkForOverwriteRisk(self.pe1Out, self.sampleName, self.clobber)
        if not self.pairedEnd:
            self.pe2Out = ""
        else:
            self.pe2Out = self.outputDirectory + self.sampleName + ".pe2.sai"
            self.clobber = runnerSupport.checkForOverwriteRisk(self.pe2Out, self.sampleName, self.clobber)
        self.samOut = self.outputDirectory + self.sampleName + ".sam"
        self.clobber = runnerSupport.checkForOverwriteRisk(self.samOut, self.sampleName, self.clobber)
                
    def makeBWAlignCommand(self, inputFile, outputFile):
        import runnerSupport
        flagValues = {"-q" : self.qualityThreshold,
                      "-l" : self.seedLength,
                      "-k" : self.maxMismatchesPerSeed,
                      "-t" : self.threads,
                      "-o" : self.maximumGapOpens,
                      "-f" : outputFile}
        bwaArgs = [programPaths["bwa"], "aln", flagValues, self.refGenomeFasta, inputFile]
        argumentFormatter = runnerSupport.ArgumentFormatter(bwaArgs)
        bwaCommand = argumentFormatter.argumentString
        return bwaCommand
        
    def makeSAMCommand(self):
        import runnerSupport
        if self.pairedEnd:
            flagValues = {"-P" : True}  #This is valid and will load everything into memory
                          # "-T" : True,  #Can't find this as a valid argument.  Copied from Catie's email.  Remove if this causes trouble.     It did.
                          # "-t" : self.threads}  #Copied from email as above.  Remove if it causes trouble.  I don't think sampe or samse can multithread operations.  If they can't, consider removing this, as it will slow our resource allocation times.
            bwaArgs = [programPaths["bwa"], "sampe", flagValues, self.refGenomeFasta, self.pe1Out, self.pe2Out, self.pe1, self.pe2, ">", self.samOut]
            argumentFormatter = runnerSupport.ArgumentFormatter(bwaArgs)
            bwaCommand = argumentFormatter.argumentString
            return bwaCommand
        else:
            # flagValues = {"-T" : True,  #Can't find this as a valid argument.  Copied from Catie's email.  Remove if this causes trouble.
            #               "-t" : self.threads}  #Copied from email as above.  Remove if it causes trouble.  I don't think sampe or samse can multithread operations.  If they can't, consider removing this, as it will slow our resource allocation times.
            bwaArgs = [programPaths["bwa"], "samse", self.refGenomeFasta, self.pe1, self.pe1Out, ">", self.samOut]
            argumentFormatter = runnerSupport.ArgumentFormatter(bwaArgs)
            bwaCommand = argumentFormatter.argumentString
            return bwaCommand

class BWAmem(object):
    
    def __init__(self, sampleName, refGenomeFasta, pairedEnd1File, pairedEnd2File = "", cores = 1, qualityThreshold = 5, seedLength = 32, maxMismatchesPerSeed = 2, maximumGapOpens = 1, clobber = None, outputDirectory = ""):
        import os
        import runnerSupport
        if not outputDirectory:
            self.outputDirectory = ""
        else:
            if not os.path.isdir(outputDirectory):
                raise RuntimeError("Output directory %s does not exist.  Please make the directory before creating jobs." %(outputDirectory))
            if not outputDirectory.endswith(os.sep):
                self.outputDirectory = outputDirectory + os.sep
            else:
                self.outputDirectory = outputDirectory
        self.sampleName = sampleName
        self.pe1 = pairedEnd1File
        self.pe2 = pairedEnd2File
        self.cores = cores
        self.threads = cores
        self.qualityThreshold = qualityThreshold
        self.seedLength = seedLength
        self.maxMismatchesPerSeed = maxMismatchesPerSeed
        self.maximumGapOpens = maximumGapOpens
        self.refGenomeFasta = refGenomeFasta
        self.clobber = clobber
        #SANITY CHECK ALL THE THINGS
        if not type(self.pe1) == str or not self.pe1:
            raise RuntimeError("Paired end 1 file must be specified and must be a string.  Value passed was \"%s\"" %(self.pe1))
        runnerSupport.checkForRequiredFile(self.pe1, "paired end 1 file")
        if not type(self.pe2) == str:
            raise RuntimeError("Paired end 2 file must be specified and must be a string.  Value passed was \"%s\"" %(self.pe2))
        if self.pe2:
            runnerSupport.checkForRequiredFile(self.pe2, "paired end 2 file")
            self.pairedEnd = True
        else:
            self.pairedEnd = False
        if not type(self.threads) == int and self.threads > 0:
            raise RuntimeError("Thread count must be a positive integer. %s was passed." %(self.threads))
        if not type(self.qualityThreshold) == int and self.qualityThreshold > 0:
            raise RuntimeError("Thread count must be a positive integer. %s was passed." %(self.qualityThreshold))
        if not type(self.seedLength) == int and self.seedLength > 0:
            raise RuntimeError("Thread count must be a positive integer. %s was passed." %(self.seedLength))
        if not type(self.maxMismatchesPerSeed) == int and self.maxMismatchesPerSeed > 0:
            raise RuntimeError("Thread count must be a positive integer. %s was passed." %(self.maxMismatchesPerSeed))
        if not type(self.maximumGapOpens) == int and self.maximumGapOpens > 0:
            raise RuntimeError("Thread count must be a positive integer. %s was passed." %(self.maximumGapOpens))
        if not type(self.refGenomeFasta) == str or not self.refGenomeFasta:
            raise RuntimeError("Reference genome fasta must be passed as a string. %s was passed." %(self.refGenomeFasta))
        self.checkRefGenome()
        self.makeAndCheckOutputFileNames()
        #DONE SANITY CHECKING ALL THE THINGS. FOR NOW.
        self.bwaCommand = self.makeBWAlignCommand()
    
    def checkRefGenome(self):
        import runnerSupport
        import os
        runnerSupport.checkForRequiredFile(self.refGenomeFasta, "reference genome file")
        runnerSupport.checkForRequiredFile(self.refGenomeFasta + ".fai", "reference genome fasta index", "Please move the index to this location or use samtools faidx to create one.")
        bwaIndexFileExtensions = [".amb", ".ann", ".bwt", ".pac", ".sa"]
        for extension in bwaIndexFileExtensions:
            indexFile = self.refGenomeFasta + extension
            runnerSupport.checkForRequiredFile(indexFile, "one or more BWA index files", "Please move the index to this location or use bwa index to (re)create one.")
    
    def makeAndCheckOutputFileNames(self):
        import runnerSupport
        self.samOut = self.outputDirectory + self.sampleName + ".sam"
        self.clobber = runnerSupport.checkForOverwriteRisk(self.samOut, self.sampleName, self.clobber)
                
    def makeBWAlignCommand(self):
        import runnerSupport
        flagValues = {}
        if self.pairedEnd:
            inputs = [self.pe1, self.pe2]
        else:
            inputs = self.pe1
        bwaArgs = [programPaths["bwa"], "mem", "-t %s" %self.cores,  self.refGenomeFasta, inputs, ">", self.samOut]
        argumentFormatter = runnerSupport.ArgumentFormatter(bwaArgs, delimiter = " ")
        bwaCommand = argumentFormatter.argumentString
        return bwaCommand

class MPileup(object):
    
    def __init__(self, sampleName, bamFile, refGenomeFasta, disablePerBaseAlignmentQuality = False, minBaseQuality = False, maxDepth = False, countOrphans = False, bgzip = False, clobber = False, outputDirectory = ""):
        import runnerSupport
        self.sampleName = sampleName
        self.bgzip = bgzip
        if not outputDirectory:
            self.outputDirectory = ""
        else:
            if not os.path.isdir(outputDirectory):
                raise RuntimeError("Output directory %s does not exist.  Please make the directory before creating jobs." %(outputDirectory))
            if not outputDirectory.endswith(os.sep):
                self.outputDirectory = outputDirectory + os.sep
            else:
                self.outputDirectory = outputDirectory
        self.bamFile = bamFile
        self.refGenomeFasta = refGenomeFasta
        self.disablePerBaseAlignmentQuality = disablePerBaseAlignmentQuality
        self.minBaseQuality = minBaseQuality
        self.maxDepth = maxDepth
        self.countOrphans = countOrphans
        self.clobber = clobber
        #SANITY CHECKING
        if not type(self.bamFile) == str:
            raise RuntimeError("BAM file name should be string. Passed %s" %(self.bamFile))
        runnerSupport.checkForRequiredFile(self.bamFile, "BAM file to analyze")
        if not type(self.refGenomeFasta) == str:
            raise RuntimeError("Reference genome fasta file name should be a string. Passed %s" %(self.refGenomeFasta))
        runnerSupport.checkForRequiredFile(self.refGenomeFasta, "reference genome fasta")
        runnerSupport.checkForRequiredFile(self.refGenomeFasta + ".fai", "reference genome fasta index file", "Please move the index file to this location or create one using samtools faidx.")
        if not type(self.minBaseQuality) == int and self.minBaseQuality > 0:
            raise RuntimeError("Minimum base quality must be a positive integer.")
        if not type(self.maxDepth) == int and self.maxDepth > 0:
            raise RuntimeError("Maximum depth must be a positive integer.")
        if not type(self.countOrphans) == bool:
            raise RuntimeError("Count Orphans must be a boolean argument.")
        #Done sanity checking
        self.makeAndCheckOutputFileNames()
        self.mPileupCommand = self.createSamtoolsCommand()
        
    def makeAndCheckOutputFileNames(self):
        import runnerSupport
        self.mPileupOut = self.outputDirectory + self.sampleName + ".mpileup"
        if self.bgzip:
            self.mPileupOut += ".bgz"
        self.clobber = runnerSupport.checkForOverwriteRisk(self.mPileupOut, self.sampleName, self.clobber)
        
    def createSamtoolsCommand(self):
        import runnerSupport
        flagValues = {"-B" : self.disablePerBaseAlignmentQuality,
                      "-Q" : self.minBaseQuality,
                      "-d" : self.maxDepth,
                      "-A" : self.countOrphans,
                      "-f" : self.refGenomeFasta}
        if self.bgzip:
            outputTee = ">" + programPaths["bgzip"] + "-c > " + self.mPileupOut
        else:
            outputTee = ">" + self.mPileupOut
        samtoolsArgs = [programPaths["samtools"], "mpileup", flagValues, self.bamFile, outputTee]
        argumentFormatter = runnerSupport.ArgumentFormatter(samtoolsArgs)
        samtoolsCommand = argumentFormatter.argumentString
        return samtoolsCommand
    
class ExtractVariantsTumor(object):
    
    def __init__(self, sampleName, pileupInput, clobber = False, minSupport = 0, requireDoubleStranded = False, outputDirectory = ""):
        import runnerSupport
        self.sampleName = sampleName
        self.requireDoubleStranded = requireDoubleStranded
        if not outputDirectory:
            self.outputDirectory = ""
        else:
            if not os.path.isdir(outputDirectory):
                raise RuntimeError("Output directory %s does not exist.  Please make the directory before creating jobs." %(outputDirectory))
            if not outputDirectory.endswith(os.sep):
                self.outputDirectory = outputDirectory + os.sep
            else:
                self.outputDirectory = outputDirectory
        self.pileupInput = pileupInput
        self.minSupport = minSupport
        self.clobber = clobber
        #sanity checking
        if not type(self.pileupInput) == str:
            raise RuntimeError("Input VCF file name must be passed as a string. Passed: %s" %(self.pileupInput))
        runnerSupport.checkForRequiredFile(self.pileupInput, "VCF input file")
        #done sanity checking
        self.makeAndCheckOutputFileNames()
        self.extractCommand = self.createPileupToVCFCommand()    
    
    def makeAndCheckOutputFileNames(self):
        import runnerSupport
        self.variantsOut = self.outputDirectory + runnerSupport.stripDirectoryAndExtension(self.pileupInput) + ".variants"
        self.targetList = self.outputDirectory + runnerSupport.stripDirectoryAndExtension(self.pileupInput) + ".targets"
        self.clobber = runnerSupport.checkForOverwriteRisk(self.variantsOut, self.sampleName, self.clobber)
        self.clobber = runnerSupport.checkForOverwriteRisk(self.targetList, self.sampleName, self.clobber)
    
    def createPileupToVCFCommand(self):
        import runnerSupport
        flagValues = {"-f" : self.pileupInput,
                      "-o" : self.variantsOut,
                      "-n" : self.minSupport,
                      "-t" : self.targetList,
                      "-d" : self.requireDoubleStranded}
        pileupCommandArgs = [programPaths["python3"], programPaths["extractVariants"], flagValues]
        argumentFormatter = runnerSupport.ArgumentFormatter(pileupCommandArgs)
        pileupCommand = argumentFormatter.argumentString
        return pileupCommand
    
class ExtractVariantsNormal(object):
    
    def __init__(self, sampleName, pileupInput, targetList, comparison, minSupport = 0, clobber = False, outputDirectory = ""):
        import runnerSupport
        self.sampleName = sampleName
        self.minSupport = minSupport
        if not outputDirectory:
            self.outputDirectory = ""
        else:
            if not os.path.isdir(outputDirectory):
                raise RuntimeError("Output directory %s does not exist.  Please make the directory before creating jobs." %(outputDirectory))
            if not outputDirectory.endswith(os.sep):
                self.outputDirectory = outputDirectory + os.sep
            else:
                self.outputDirectory = outputDirectory
        if not comparison:
            self.comparison = ""
        else:
            self.comparison = "." + comparison
        self.pileupInput = pileupInput
        self.targetList = targetList
        self.clobber = clobber
        #sanity checking
        if not type(self.pileupInput) == str:
            raise RuntimeError("Input VCF file name must be passed as a string. Passed: %s" %(self.pileupInput))
        runnerSupport.checkForRequiredFile(self.pileupInput, "VCF input file")
        runnerSupport.checkForRequiredFile(self.targetList, "List of targets from tumor to do ")
        #done sanity checking
        self.makeAndCheckOutputFileNames()
        self.extractCommand = self.createPileupToVCFCommand()
    
    def makeAndCheckOutputFileNames(self):
        import runnerSupport
        self.variantsOut = self.outputDirectory + runnerSupport.stripDirectoryAndExtension(self.pileupInput) + self.comparison + ".variants"
        self.clobber = runnerSupport.checkForOverwriteRisk(self.variantsOut, self.sampleName, self.clobber)
    
    def createPileupToVCFCommand(self):
        import runnerSupport
        flagValues = {"-f" : self.pileupInput,
                      "-o" : self.variantsOut,
                      "-n" : self.minSupport,
                      "-m" : self.targetList}
        pileupCommandArgs = [programPaths["python3"], programPaths["extractVariants"], flagValues]
        argumentFormatter = runnerSupport.ArgumentFormatter(pileupCommandArgs)
        pileupCommand = argumentFormatter.argumentString
        return pileupCommand
    
class CombineExtractedVariants(object):
    
    def __init__(self, sampleName, tumorFileName, normalFileName, comparison, clobber = False, outputDirectory = ""):
        import runnerSupport
        self.sampleName = sampleName
        if not outputDirectory:
            self.outputDirectory = ""
        else:
            if not os.path.isdir(outputDirectory):
                raise RuntimeError("Output directory %s does not exist.  Please make the directory before creating jobs." %(outputDirectory))
            if not outputDirectory.endswith(os.sep):
                self.outputDirectory = outputDirectory + os.sep
            else:
                self.outputDirectory = outputDirectory
        if not comparison:
            self.comparison = ""
        else:
            self.comparison = "." + comparison
        self.tumorFile = tumorFileName
        self.normalFile = normalFileName
        self.clobber = clobber
        #sanity checking
        runnerSupport.checkForRequiredFile(self.tumorFile, "Tumor input file")
        runnerSupport.checkForRequiredFile(self.normalFile, "Normal input file")
        #done sanity checking
        self.makeAndCheckOutputFileNames()
        self.combineCommand = self.createPileupToVCFCommand()    
    
    def makeAndCheckOutputFileNames(self):
        import runnerSupport
        self.vcfOut = self.outputDirectory + runnerSupport.stripDirectoryAndExtension(self.normalFile) + self.comparison + ".vcf"
        self.clobber = runnerSupport.checkForOverwriteRisk(self.vcfOut, self.sampleName, self.clobber)
    
    def createPileupToVCFCommand(self):
        import runnerSupport
        flagValues = {"-t" : self.tumorFile,
                      "-n" : self.normalFile,
                      "-o" : self.vcfOut}
        combineCommandArgs = [programPaths["python3"], programPaths["combineVariants"], flagValues]
        argumentFormatter = runnerSupport.ArgumentFormatter(combineCommandArgs)
        combineCommand = argumentFormatter.argumentString
        return combineCommand

class Tabix(object):
    
    def __init__(self, sampleName, pileupIn, sequenceColumn = 1, beginColumn = 2, endColumn = 2, clobber = False, outputDirectory = ""):
        import runnerSupport
        self.sampleName = sampleName
        self.pileupIn = pileupIn
        self.sequenceColumn = sequenceColumn
        self.beginColumn = beginColumn
        self.endColumn = endColumn
        #SANITY TEST ALL THE THINGS
        runnerSupport.checkForRequiredFile(self.bamIn, "BAM file for depth check")
        if not type(self.sequenceColumn) == int:
            raise RuntimeError("Sequence column must be an integer.")
        if not type(self.beginColumn) == int:
            raise RuntimeError("Begin column must be an integer.")
        if not type(self.endColumn) == int:
            raise RuntimeError("End column must be an integer.")
        runnerSupport.checkForRequiredFile(self.pileupIn, "Pileup file to index")
        #DONE SANITY CHECKING. FOR NOW.
        self.makeAndCheckOutputFileNames()
        self.tabixCommand = self.createTabixCommand()
        
    def makeAndCheckOutputFileNames(self):
        import runnerSupport
        self.indexOut = self.outputDirectory + runnerSupport.stripDirectoryAndExtension(self.pileupIn) + ".bgz.tbi" 
        self.clobber = runnerSupport.checkForOverwriteRisk(self.indexOut, self.sampleName, self.clobber)
        
    def createTabixCommand(self):
        import runnerSupport
        flagValues = {"-s" : self.sequenceColumn,
                      "-b" : self.beginColumn,
                      "-e" : self.endColumn}
        gatkArgs = [programPaths["tabix"], flagValues, self.pileupIn]
        argumentFormatter = runnerSupport.ArgumentFormatter(gatkArgs)
        depthCommand = argumentFormatter.argumentString
        return (depthCommand)
    
class ViewSAMtoBAM(object):
    
    def __init__(self, sampleName, samFile, refGenomeFasta, samInput = True, bamOutput = True, includeHeaders = True, clobber = False, outputDirectory = ""):
        import runnerSupport
        self.samInput = samInput
        self.bamOutput = bamOutput
        self.includeHeaders = includeHeaders
        self.sampleName = sampleName
        if not outputDirectory:
            self.outputDirectory = ""
        else:
            if not os.path.isdir(outputDirectory):
                raise RuntimeError("Output directory %s does not exist.  Please make the directory before creating jobs." %(outputDirectory))
            if not outputDirectory.endswith(os.sep):
                self.outputDirectory = outputDirectory + os.sep
            else:
                self.outputDirectory = outputDirectory
        self.samFile = samFile
        self.refGenomeFasta = refGenomeFasta
        self.clobber = clobber
        #SANITY CHECKING
        runnerSupport.checkTypes([self.samInput, self.bamOutput, self.includeHeaders], bool)
        if not type(self.samFile) == str:
            raise RuntimeError("BAM file name should be string. Passed %s" %(self.samFile))
        runnerSupport.checkForRequiredFile(self.samFile, "SAM file to convert")
        if not type(self.refGenomeFasta) == str:
            raise RuntimeError("Reference genome fasta file name should be a string. Passed %s" %(self.refGenomeFasta))
        runnerSupport.checkForRequiredFile(self.refGenomeFasta, "reference genome fasta")
        runnerSupport.checkForRequiredFile(self.refGenomeFasta + ".fai", "reference genome fasta index file", "Please move the index file to this location or create one using samtools faidx.")
        #Done sanity checking
        self.makeAndCheckOutputFileNames()
        self.viewCommand = self.createSamtoolsCommand()
        
    def makeAndCheckOutputFileNames(self):
        import runnerSupport
        self.bamOut = self.outputDirectory + self.sampleName + ".bam"
        self.clobber = runnerSupport.checkForOverwriteRisk(self.bamOut, self.sampleName, self.clobber)
        
    def createSamtoolsCommand(self):
        import runnerSupport
        flagValues = {"-S" : self.samInput,
                      "-b" : self.bamOutput,
                      "-h" : self.includeHeaders,
                      "-T" : self.refGenomeFasta,
                      "-o" : self.bamOut}
        samtoolsArgs = [programPaths["samtools"], "view", flagValues, self.samFile]
        argumentFormatter = runnerSupport.ArgumentFormatter(samtoolsArgs)
        samtoolsCommand = argumentFormatter.argumentString
        return samtoolsCommand
    
class GrepFilter(object):
    
    def __init__(self, sampleName, inputFile, outputFile, filterString, filterOut = False, clobber = False, outputDirectory = ""):
        import runnerSupport
        self.clobber = clobber
        self.sampleName = sampleName
        self.inputFile = inputFile
        self.outputFile = outputFile
        self.filterString = filterString
        self.filterOut = filterOut
        if not outputDirectory:
            self.outputDirectory = ""
        else:
            if not os.path.isdir(outputDirectory):
                raise RuntimeError("Output directory %s does not exist.  Please make the directory before creating jobs." %(outputDirectory))
            if not outputDirectory.endswith(os.sep):
                self.outputDirectory = outputDirectory + os.sep
            else:
                self.outputDirectory = outputDirectory
        #SANITY CHECKING
        runnerSupport.checkTypes([self.filterOut], bool)
        runnerSupport.checkForRequiredFile(self.inputFile, "File to filter")
        #Done sanity checking
        self.makeAndCheckOutputFileNames()
        self.grepCommand = self.createGrepCommand()
        
    def makeAndCheckOutputFileNames(self):
        import runnerSupport
        self.clobber = runnerSupport.checkForOverwriteRisk(self.outputFile, self.sampleName, self.clobber)
        
    def createGrepCommand(self):
        import runnerSupport
        flagValues = {"-v" : self.filterOut}
        grepArgs = ["grep", flagValues, self.filterString, self.inputFile, ">", self.outputFile]
        argumentFormatter = runnerSupport.ArgumentFormatter(grepArgs)
        grepCommand = argumentFormatter.argumentString
        return grepCommand
    
class Varscan(object):
    
    def __init__(self, sampleName, normalPileup, tumorPileup, minimumTumorVariantFrequency = 0.05, tumorPurity = 1.0, homozygousFrequency = 0.75, clobber = False, outputDirectory = ""):
        import runnerSupport
        self.clobber = clobber
        self.sampleName = sampleName
        self.normalPileup = normalPileup
        self.tumorPileup = tumorPileup
        self.minimumTumorVariantFrequency = minimumTumorVariantFrequency
        self.tumorPurity = tumorPurity
        self.homozygousFrequency = homozygousFrequency
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
        runnerSupport.checkForRequiredFile(self.normalPileup, "Normal mPileup file")
        runnerSupport.checkForRequiredFile(self.tumorPileup, "Tumor mPileup file")
        runnerSupport.checkTypes([self.minimumTumorVariantFrequency, self.tumorPurity, self.homozygousFrequency], float)
        #DONE SANITY CHECKING. FOR NOW.
        self.makeAndCheckOutputFileNames()
        self.varscanSomaticCommand, self.varscanProcessSNPCommand, self.varscanProcessIndelCommand = self.createVarscanCommand()
        
    def makeAndCheckOutputFileNames(self):
        import runnerSupport
        self.baseName = self.outputDirectory + self.sampleName + ".varscan"
        self.snpOut = self.outputDirectory + self.sampleName + ".varscan.snp"
        self.indelOut = self.outputDirectory + self.sampleName + ".varscan.indel"
        self.hcSomaticSNPOut = self.snpOut + ".Somatic.hc"
        self.hcSomaticIndelOut = self.indelOut + ".Somatic.hc"
        self.clobber = runnerSupport.checkForOverwriteRisk(self.snpOut, self.sampleName, self.clobber)
        self.clobber = runnerSupport.checkForOverwriteRisk(self.indelOut, self.sampleName, self.clobber)
        self.clobber = runnerSupport.checkForOverwriteRisk(self.hcSomaticSNPOut, self.sampleName, self.clobber)
        self.clobber = runnerSupport.checkForOverwriteRisk(self.hcSomaticIndelOut, self.sampleName, self.clobber)

    def createVarscanCommand(self):
        import runnerSupport
        flagValuesSomatic = {"--min-var-freq" : self.minimumTumorVariantFrequency,
                             "--tumor-purity" : self.tumorPurity,
                             "--min-freq-for-hom" : self.homozygousFrequency}
        flagValuesProcess = {"--min-tumor-freq" : self.minimumTumorVariantFrequency}
        somaticArgs = [programPaths["java"], "-Xmx1g", "-jar", programPaths["varscan"], "somatic", self.normalPileup, self.tumorPileup, self.baseName , flagValuesSomatic]
        snpArgs = [programPaths["java"], "-Xmx1g", "-jar", programPaths["varscan"], "processSomatic", self.snpOut, flagValuesProcess]
        indelArgs = [programPaths["java"], "-Xmx1g", "-jar", programPaths["varscan"], "processSomatic", self.indelOut, flagValuesProcess]
        argumentFormatter = runnerSupport.ArgumentFormatter(somaticArgs)
        somaticCommand = argumentFormatter.argumentString
        argumentFormatter = runnerSupport.ArgumentFormatter(snpArgs)
        snpCommand = argumentFormatter.argumentString
        argumentFormatter = runnerSupport.ArgumentFormatter(indelArgs)
        indelCommand = argumentFormatter.argumentString
        return (somaticCommand, snpCommand, indelCommand)
    
class BAMReadCount(object):
    
    def __init__(self, sampleName, bamIn, inputPositions, refGenomeFasta, maxWarnings = 10, clobber = False, outputDirectory = ""):
        import runnerSupport
        self.bamIn = bamIn
        self.clobber = clobber
        self.refGenomeFasta = refGenomeFasta
        self.sampleName = sampleName
        self.maxWarnings = maxWarnings
        self.inputPositions = inputPositions
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
        runnerSupport.checkForRequiredFile(self.inputPositions, "Input position list")
        runnerSupport.checkTypes([self.maxWarnings], int)
        runnerSupport.checkForRequiredFile(self.bamIn, "BAM file to analyze for read counts")
        #DONE SANITY CHECKING. FOR NOW.
        self.checkRefGenome()
        self.makeAndCheckOutputFileNames()
        self.bamReadcountCommand = self.createBAMReadcountCommand()
        
    def makeAndCheckOutputFileNames(self):
        import runnerSupport
        self.readcountOut = self.outputDirectory + runnerSupport.stripDirectory(self.inputPositions) + ".readcount"
        self.clobber = runnerSupport.checkForOverwriteRisk(self.readcountOut, self.sampleName, self.clobber)

    def checkRefGenome(self):
        import runnerSupport
        import os
        runnerSupport.checkForRequiredFile(self.refGenomeFasta, "reference genome file")
        runnerSupport.checkForRequiredFile(self.refGenomeFasta + ".fai", "reference genome fasta index", "Please move the index to this location or use samtools faidx to create one.")

    def createBAMReadcountCommand(self):
        import runnerSupport
        flagValues = {"-w" : self.maxWarnings,
                      "-f" : self.refGenomeFasta,
                      "-l" : self.inputPositions}
        args = [programPaths["bam-readcount"], flagValues, self.bamIn, ">" , self.readcountOut]
        argumentFormatter = runnerSupport.ArgumentFormatter(args)
        bamReadcountCommand = argumentFormatter.argumentString
        return bamReadcountCommand
    
class VarScanPositions(object):
    
    def __init__(self, sampleName, varScanData, clobber = False, outputDirectory = ""):
        import runnerSupport
        self.varScanData = varScanData
        self.clobber = clobber
        self.sampleName = sampleName
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
        runnerSupport.checkForRequiredFile(self.varScanData, "VarScan output data")
        #DONE SANITY CHECKING. FOR NOW.
        self.makeAndCheckOutputFileNames()
        self.positionPullerCommand = self.createPositionPullerCommand()
        
    def makeAndCheckOutputFileNames(self):
        import runnerSupport
        self.positionsOut = self.outputDirectory + runnerSupport.stripDirectory(self.varScanData) + ".positions"
        self.clobber = runnerSupport.checkForOverwriteRisk(self.positionsOut, self.sampleName, self.clobber)

    def createPositionPullerCommand(self):
        import runnerSupport
        flagValues = {"-f" : self.varScanData,
                      "-o" : self.positionsOut}
        args = [programPaths["python3"], programPaths["varScanPositions"], flagValues]
        argumentFormatter = runnerSupport.ArgumentFormatter(args)
        positionPullerCommand = argumentFormatter.argumentString
        return positionPullerCommand
    
class VCFReader(object):
    
    def __init__(self, sampleName, vcfData, tumorName, normalName, maxPValue = 0.05, minDepth = 10, clobber = False, outputDirectory = ""):
        import runnerSupport
        self.vcfData = vcfData
        self.tumorName = tumorName
        self.normalName = normalName
        self.clobber = clobber
        self.maxPValue = maxPValue
        self.minDepth = minDepth
        self.sampleName = sampleName
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
        runnerSupport.checkForRequiredFile(self.vcfData, "VCF (variant call file)")
        runnerSupport.checkTypes([self.minDepth], [int, bool])
        runnerSupport.checkTypes([self.maxPValue], [float, int, bool])
        #DONE SANITY CHECKING. FOR NOW.
        self.makeAndCheckOutputFileNames()
        self.vcfReaderCommand = self.createVCFReaderCommand()
        
    def makeAndCheckOutputFileNames(self):
        import runnerSupport
        self.acceptedOut = self.outputDirectory + runnerSupport.stripDirectory(self.vcfData) + ".haplotypeCaller.accepted.pkl"
        self.clobber = runnerSupport.checkForOverwriteRisk(self.acceptedOut, self.sampleName, self.clobber)

    def createVCFReaderCommand(self):
        import runnerSupport
        flagValues = {"-f" : self.vcfData,
                      "-o" : self.acceptedOut,
                      "-t" : self.tumorName,
                      "-n" : self.normalName,
                      "-p" : self.maxPValue,
                      "-d" : self.minDepth}
        args = [programPaths["python3"], programPaths["vcfReader"], flagValues]
        argumentFormatter = runnerSupport.ArgumentFormatter(args)
        vcfReaderCommand = argumentFormatter.argumentString
        return vcfReaderCommand
    
class MutectReader(object):
    
    def __init__(self, sampleName, mutectData, minDepth = 10, clobber = False, outputDirectory = ""):
        import runnerSupport
        self.mutectData = mutectData
        self.clobber = clobber
        self.minDepth = minDepth
        self.sampleName = sampleName
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
        runnerSupport.checkForRequiredFile(self.mutectData, "Mutect output data")
        runnerSupport.checkTypes([self.minDepth], [int, bool])
        #DONE SANITY CHECKING. FOR NOW.
        self.makeAndCheckOutputFileNames()
        self.mutectReaderCommand = self.createMutectReaderCommand()
        
    def makeAndCheckOutputFileNames(self):
        import runnerSupport
        self.acceptedOut = self.outputDirectory + runnerSupport.stripDirectory(self.mutectData) + ".mutect1.accepted.pkl"
        self.clobber = runnerSupport.checkForOverwriteRisk(self.acceptedOut, self.sampleName, self.clobber)

    def createMutectReaderCommand(self):
        import runnerSupport
        flagValues = {"-f" : self.mutectData,
                      "-o" : self.acceptedOut,
                      "-d" : self.minDepth}
        args = [programPaths["python3"], programPaths["mutectReader"], flagValues]
        argumentFormatter = runnerSupport.ArgumentFormatter(args)
        mutectReaderCommand = argumentFormatter.argumentString
        return mutectReaderCommand
    
class VarScanReader(object):
    
    def __init__(self, sampleName, varScanData, minDepth = 10, clobber = False, outputDirectory = ""):
        import runnerSupport
        self.varScanData = varScanData
        self.clobber = clobber
        self.minDepth = minDepth
        self.sampleName = sampleName
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
        runnerSupport.checkForRequiredFile(self.varScanData, "VarScan output data")
        runnerSupport.checkTypes([self.minDepth], [int, bool])
        #DONE SANITY CHECKING. FOR NOW.
        self.makeAndCheckOutputFileNames()
        self.varScanReaderCommand = self.createVarScanReaderCommand()
        
    def makeAndCheckOutputFileNames(self):
        import runnerSupport
        self.acceptedOut = self.outputDirectory + runnerSupport.stripDirectory(self.varScanData) + ".varScan.accepted.pkl"
        self.clobber = runnerSupport.checkForOverwriteRisk(self.acceptedOut, self.sampleName, self.clobber)

    def createVarScanReaderCommand(self):
        import runnerSupport
        flagValues = {"-f" : self.varScanData,
                      "-o" : self.acceptedOut,
                      "-d" : self.minDepth}
        args = [programPaths["python3"], programPaths["varScanReader"], flagValues]
        argumentFormatter = runnerSupport.ArgumentFormatter(args)
        varScanReaderCommand = argumentFormatter.argumentString
        return varScanReaderCommand

class VarScanFPFilter(object):
    
    def __init__(self, sampleName, varScanData, bamReadCounts, clobber = False, outputDirectory = ""):
        import runnerSupport
        self.varScanData = varScanData
        self.bamReadCounts = bamReadCounts
        self.clobber = clobber
        self.sampleName = sampleName
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
        runnerSupport.checkForRequiredFile(self.varScanData, "VarScan output data")
        runnerSupport.checkForRequiredFile(self.bamReadCounts, "BAM read count data")
        #DONE SANITY CHECKING. FOR NOW.
        self.makeAndCheckOutputFileNames()
        self.fpFilterCommand = self.createFPFilterCommand()
        
    def makeAndCheckOutputFileNames(self):
        import runnerSupport
        self.baseFileName = self.outputDirectory + runnerSupport.stripDirectory(self.varScanData)
        self.passOut = self.baseFileName + ".pass"
        self.clobber = runnerSupport.checkForOverwriteRisk(self.passOut, self.sampleName, self.clobber)

    def createFPFilterCommand(self):
        import runnerSupport
        flagValues = {"--output-basename" : self.baseFileName}
        args = [programPaths["perl"], programPaths["varScanFPFilter"], self.varScanData, self.bamReadCounts, flagValues]
        argumentFormatter = runnerSupport.ArgumentFormatter(args)
        fpFilterCommand = argumentFormatter.argumentString
        return fpFilterCommand
    
class VariantCombine(object):
    
    def __init__(self, sampleName, variantDataList, minHits = False, maxHits = False, clobber = False, outputDirectory = ""):
        import runnerSupport
        if type(variantDataList) == str:
            variantDataList = [variantDataList]
        self.variantDataList = variantDataList
        self.minHits = minHits
        self.maxHits = maxHits
        self.clobber = clobber
        self.sampleName = sampleName
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
        for file in variantDataList:
            runnerSupport.checkForRequiredFile(file.split(",")[-1], "Variant data file")
        runnerSupport.checkTypes([minHits, maxHits], (bool, int))
        #DONE SANITY CHECKING. FOR NOW.
        self.makeAndCheckOutputFileNames()
        self.variantCombineCommand = self.createCombineCommand()
        
    def makeAndCheckOutputFileNames(self):
        import runnerSupport
        self.acceptedSomaticVariants = self.outputDirectory + self.sampleName + ".acceptedSomatics.txt"
        self.clobber = runnerSupport.checkForOverwriteRisk(self.acceptedSomaticVariants, self.sampleName, self.clobber)

    def createCombineCommand(self):
        import runnerSupport
        flagValues = {"flaggedlist" : ["-v", self.variantDataList],
                      "-x" : self.maxHits,
                      "-m" : self.minHits,
                      "-o" : self.acceptedSomaticVariants}
        args = [programPaths["python3"], programPaths["variantCombine"], flagValues]
        argumentFormatter = runnerSupport.ArgumentFormatter(args)
        variantCombineCommand = argumentFormatter.argumentString
        return variantCombineCommand
    
class Hisat2Align(object):
    
    def __init__(self, sampleName, refGenomePrefix, pairedEnd1File, pairedEnd2File = "", cores = 1, clobber = None, outputDirectory = ""):
        import os
        import runnerSupport
        if not outputDirectory:
            self.outputDirectory = ""
        else:
            if not os.path.isdir(outputDirectory):
                raise RuntimeError("Output directory %s does not exist.  Please make the directory before creating jobs." %(outputDirectory))
            if not outputDirectory.endswith(os.sep):
                self.outputDirectory = outputDirectory + os.sep
            else:
                self.outputDirectory = outputDirectory
        self.sampleName = sampleName
        self.pe1 = pairedEnd1File
        self.pe2 = pairedEnd2File
        self.cores = cores
        self.threads = cores
        self.refGenomePrefix = refGenomePrefix
        self.clobber = clobber
        #SANITY CHECK ALL THE THINGS
        if not type(self.pe1) == str or not self.pe1:
            raise RuntimeError("Paired end 1 file must be specified and must be a string.  Value passed was \"%s\"" %(self.pe1))
        runnerSupport.checkForRequiredFile(self.pe1, "paired end 1 file")
        if not type(self.pe2) == str:
            raise RuntimeError("Paired end 2 file must be specified and must be a string.  Value passed was \"%s\"" %(self.pe2))
        if self.pe2:
            runnerSupport.checkForRequiredFile(self.pe2, "paired end 2 file")
            self.pairedEnd = True
        else:
            self.pairedEnd = False
        if not type(self.threads) == int and self.threads > 0:
            raise RuntimeError("Thread count must be a positive integer. %s was passed." %(self.threads))
        if not type(self.refGenomePrefix) == str or not self.refGenomePrefix:
            raise RuntimeError("Reference genome fasta must be passed as a string. %s was passed." %(self.refGenomePrefix))
        self.checkRefGenome()
        self.makeAndCheckOutputFileNames()
        #DONE SANITY CHECKING ALL THE THINGS. FOR NOW.
        self.hisat2Command = self.makeHisat2Command()
    
    def checkRefGenome(self):
        import runnerSupport
        import os
        runnerSupport.checkForRequiredFile(self.refGenomePrefix + ".1.ht2", "reference genome index")
    
    def makeAndCheckOutputFileNames(self):
        import runnerSupport
        self.samOut = self.outputDirectory + self.sampleName + ".sam"
        self.clobber = runnerSupport.checkForOverwriteRisk(self.samOut, self.sampleName, self.clobber)
                
    def makeHisat2Command(self):
        import runnerSupport
        flagValues = {"-S" : self.samOut,
                      "-1" : self.pe1,
                      "-2" : self.pe2,
                      "-p" : self.cores}
        hisat2Args = [programPaths["hisat2"], self.refGenomePrefix, flagValues]
        argumentFormatter = runnerSupport.ArgumentFormatter(hisat2Args)
        hisat2Command = argumentFormatter.argumentString
        return hisat2Command
    
class GetRNASupport(object):
    
    def __init__(self, sampleName, rnaSampleName, somaticVariantPickle, rnaVariantVCF, minDifference = 10, clobber = False, outputDirectory = ""):
        import runnerSupport
        self.somaticVariantPickle = somaticVariantPickle
        self.rnaVariantVCF = rnaVariantVCF
        self.minDifference = 10
        self.clobber = clobber
        self.sampleName = sampleName
        self.rnaSampleName = rnaSampleName
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
        runnerSupport.checkForRequiredFile(somaticVariantPickle, "Pickle of filtered somatic variants")
        runnerSupport.checkForRequiredFile(rnaVariantVCF, "VCF of observed RNA variants")
        runnerSupport.checkTypes(minDifference, (bool, int))
        #DONE SANITY CHECKING. FOR NOW.
        self.makeAndCheckOutputFileNames()
        self.getRNASupportCommand = self.createRNASupportCommand()
        
    def makeAndCheckOutputFileNames(self):
        import runnerSupport
        self.rnaSupportOut = self.outputDirectory + runnerSupport.stripDirectoryAndExtension(self.somaticVariantPickle) + ".rnaSupport.txt"
        self.clobber = runnerSupport.checkForOverwriteRisk(self.rnaSupportOut, self.sampleName, self.clobber)

    def createRNASupportCommand(self):
        import runnerSupport
        flagValues = {"-v" : self.rnaVariantVCF,
                      "-s" : self.somaticVariantPickle,
                      "-r" : self.rnaSampleName,
                      "-m" : self.minDifference,
                      "-o" : self.rnaSupportOut}
        args = [programPaths["python3"], programPaths["getRNASupport"], flagValues]
        argumentFormatter = runnerSupport.ArgumentFormatter(args)
        getRNASupportCommand = argumentFormatter.argumentString
        return getRNASupportCommand