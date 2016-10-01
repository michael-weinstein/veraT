#!/usr/bin/env python3

import os
global programPaths
programPaths = {"bwa" : os.getcwd() + "/bin/bwa-0.7.15/bwa",
                "java" : os.getcwd() + "/bin/jre1.8.0_77/bin/java",
                "picard" : os.getcwd() + "/bin/picard-tools-2.1.1/picard.jar",
                "samtools" : "/u/local/apps/samtools/1.2/gcc-4.4.7/bin/samtools",
                "extractVariants" : os.getcwd() + "/analysisScripts/extractVariants.py",
                "combineVariants" : os.getcwd() + "/analysisScripts/combineExtractedVariants.py",
                "python3" : "/u/local/apps/python/3.4.3/bin/python3"}

def yesAnswer(question):  #asks the question passed in and returns True if the answer is yes, False if the answer is no, and keeps the user in a loop until one of those is given.  Also useful for walking students through basic logical python functions
    answer = False  #initializes the answer variable to false.  Not absolutely necessary, since it should be undefined at this point and test to false, but explicit is always better than implicit
    while not answer:  #enters the loop and stays in it until answer is equal to True
        print (question + ' (Y/N)')  #Asks the question contained in the argument passed into this subroutine
        answer = input('>>') #sets answer equal to some value input by the user
        if str(answer) == 'y' or str(answer) == 'Y' or str(answer).upper() == 'YES':  #checks if the answer is a valid yes answer
            return True  #sends back a value of True because of the yes answer
        elif str(answer) == 'n' or str(answer) == 'N' or str(answer).upper() == 'NO': #checks to see if the answer is a valid form of no
            return False  #sends back a value of False because it was not a yes answer
        else: #if the answer is not a value indicating a yes or no
            print ('Invalid response.')
            answer = False #set ansewr to false so the loop will continue until a satisfactory answer is given

def stripDirectoryAndExtension(fileName):
    import os
    pathSplit = fileName.split(os.sep)
    fileName = pathSplit[-1]
    splitName = fileName.split(".")
    splitNameNoExtension = splitName[:-1]
    nameNoExtension = ".".join(splitNameNoExtension)
    return nameNoExtension

def checkForOverwriteRisk(file, sample = "NoneGiven", clobber = False):  #This will check if the user wants to delete the original file for overwriting.  We will delete here to avoid future collisions.  This will also return their desired change to the clobber value as a boolean.
    import os
    import time
    global rerunOfPrevious
    if os.path.isfile(file):
        print("Potential collision found for file %s" %(file))
        fileHandle = open(file, 'rb')  #reading as binary because we don't know what will actually be in there for sure
        fileSample = fileHandle.read(2000)
        fileHandle.close()
        placeholderFound = False
        if len(fileSample) >= 24:
            if fileSample[0:24] == b"PLACEHOLDER FOR SAMPLE: ":
                placeholderFound = True
                originalHolder = fileSample[24:].decode()
        if placeholderFound:
            try:
                rerun = rerunOfPrevious
            except NameError:
                rerun = False
            if not rerun:
                if originalHolder == sample:
                    print("Original placeholder belonged to %s and current run belongs to the same." %(sample))
                    print("If this is not a rerun of the same sample, please do not continue, as it will likely cause terrible file collisions.")
                    if not yesAnswer("Is this a rerun?"):
                        raise RuntimeError("Stopped to avoid collisions")
                    else:
                        rerunOfPrevious = True
                        clobber = True
        if not clobber and not placeholderFound:
            if not yesAnswer(file + " already exists. Do you wish to overwrite this file and continue?"):
                raise RuntimeError("Quit to avoid overwrite of %s." %(file))
            print("Overwriting %s" %(file))
        try:
            if not clobber:
                time.sleep(7)
        except KeyboardInterrupt:
            raise KeyboardInterrupt("Deletion aborted due to keyboard interrupt")
        os.remove(file)
        if type(clobber) == type(None):
            if yesAnswer("Continue to check file collisions?"):
                clobber = False
            else:
                clobber = True
    #Writing placeholder files to monitor for collisions with future files
    touchFile = open(file, 'w')
    touchFile.write("PLACEHOLDER FOR SAMPLE: " + sample)
    touchFile.close()
    return clobber

def checkForRequiredFile(fileName, fileDescription, instructions = ""):
    import os
    if not os.path.isfile(fileName):
        if instructions:
            if not instructions.startswith(" "):
                instructions = " " + instructions
        raise FileNotFoundError("Unable to find %s.%s Expected file: %s" %(fileDescription, instructions, fileName))

class BWAlign(object):
    
    def __init__(self, sampleName, refGenomeFasta, pairedEnd1File, pairedEnd2File = "", cores = 1, qualityThreshold = 5, seedLength = 32, maxMismatchesPerSeed = 2, maximumGapOpens = 1, clobber = None, outputDirectory = ""):
        import os
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
        checkForRequiredFile(self.pe1, "paired end 1 file")
        if not type(self.pe2) == str:
            raise RuntimeError("Paired end 2 file must be specified and must be a string.  Value passed was \"%s\"" %(self.pe2))
        if self.pe2:
            checkForRequiredFile(self.pe2, "paired end 2 file")
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
        checkForRequiredFile(self.refGenomeFasta, "reference genome file")
        checkForRequiredFile(self.refGenomeFasta + ".fai", "reference genome fasta index", "Please move the index to this location or use samtools faidx to create one.")
        bwaIndexFileExtensions = [".amb", ".ann", ".bwt", ".pac", ".sa"]
        for extension in bwaIndexFileExtensions:
            indexFile = self.refGenomeFasta + extension
            checkForRequiredFile(indexFile, "one or more BWA index files", "Please move the index to this location or use bwa index to (re)create one.")
    
    def makeAndCheckOutputFileNames(self):
        self.pe1Out = self.outputDirectory + self.sampleName + ".pe1.sai"
        self.clobber = checkForOverwriteRisk(self.pe1Out, self.sampleName, self.clobber)
        if not self.pairedEnd:
            self.pe2Out = ""
        else:
            self.pe2Out = self.outputDirectory + self.sampleName + ".pe2.sai"
            self.clobber = checkForOverwriteRisk(self.pe2Out, self.sampleName, self.clobber)
        self.samOut = self.outputDirectory + self.sampleName + ".sam"
        self.clobber = checkForOverwriteRisk(self.samOut, self.sampleName, self.clobber)
                
    def makeBWAlignCommand(self, inputFile, outputFile):
        import genericRunners
        flagValues = {"-q" : self.qualityThreshold,
                      "-l" : self.seedLength,
                      "-k" : self.maxMismatchesPerSeed,
                      "-t" : self.threads,
                      "-o" : self.maximumGapOpens,
                      "-f" : outputFile}
        bwaArgs = [programPaths["bwa"], "aln", flagValues, self.refGenomeFasta, inputFile]
        argumentFormatter = genericRunners.ArgumentFormatter(bwaArgs)
        bwaCommand = argumentFormatter.argumentString
        return bwaCommand
        
    def makeSAMCommand(self):
        import genericRunners
        if self.pairedEnd:
            flagValues = {"-P" : True}  #This is valid and will load everything into memory
                          # "-T" : True,  #Can't find this as a valid argument.  Copied from Catie's email.  Remove if this causes trouble.     It did.
                          # "-t" : self.threads}  #Copied from email as above.  Remove if it causes trouble.  I don't think sampe or samse can multithread operations.  If they can't, consider removing this, as it will slow our resource allocation times.
            bwaArgs = [programPaths["bwa"], "sampe", flagValues, self.refGenomeFasta, self.pe1Out, self.pe2Out, self.pe1, self.pe2, ">", self.samOut]
            argumentFormatter = genericRunners.ArgumentFormatter(bwaArgs)
            bwaCommand = argumentFormatter.argumentString
            return bwaCommand
        else:
            # flagValues = {"-T" : True,  #Can't find this as a valid argument.  Copied from Catie's email.  Remove if this causes trouble.
            #               "-t" : self.threads}  #Copied from email as above.  Remove if it causes trouble.  I don't think sampe or samse can multithread operations.  If they can't, consider removing this, as it will slow our resource allocation times.
            bwaArgs = [programPaths["bwa"], "samse", self.refGenomeFasta, self.pe1, self.pe1Out, ">", self.samOut]
            argumentFormatter = genericRunners.ArgumentFormatter(bwaArgs)
            bwaCommand = argumentFormatter.argumentString
            return bwaCommand

class BWAmem(object):
    
    def __init__(self, sampleName, refGenomeFasta, pairedEnd1File, pairedEnd2File = "", cores = 1, qualityThreshold = 5, seedLength = 32, maxMismatchesPerSeed = 2, maximumGapOpens = 1, clobber = None, outputDirectory = ""):
        import os
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
        checkForRequiredFile(self.pe1, "paired end 1 file")
        if not type(self.pe2) == str:
            raise RuntimeError("Paired end 2 file must be specified and must be a string.  Value passed was \"%s\"" %(self.pe2))
        if self.pe2:
            checkForRequiredFile(self.pe2, "paired end 2 file")
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
        self.bwaCommand = self.makeBWACommand()
    
    def checkRefGenome(self):
        import os
        checkForRequiredFile(self.refGenomeFasta, "reference genome file")
        checkForRequiredFile(self.refGenomeFasta + ".fai", "reference genome fasta index", "Please move the index to this location or use samtools faidx to create one.")
        bwaIndexFileExtensions = [".amb", ".ann", ".bwt", ".pac", ".sa"]
        for extension in bwaIndexFileExtensions:
            indexFile = self.refGenomeFasta + extension
            checkForRequiredFile(indexFile, "one or more BWA index files", "Please move the index to this location or use bwa index to (re)create one.")
    
    def makeAndCheckOutputFileNames(self):
        self.samOut = self.outputDirectory + self.sampleName + ".sam"
        self.clobber = checkForOverwriteRisk(self.samOut, self.sampleName, self.clobber)
                
    def makeBWACommand(self):
        import genericRunners
        flagValues = {}
        if self.pairedEnd:
            inputs = [self.pe1, self.pe2]
        else:
            inputs = self.pe1
        bwaArgs = [programPaths["bwa"], "mem", self.refGenomeFasta, inputs, ">", self.samOut]
        argumentFormatter = genericRunners.ArgumentFormatter(bwaArgs, delimiter = " ")
        bwaCommand = argumentFormatter.argumentString
        return bwaCommand

class SAMtoBAM(object):
    
    def __init__(self, sampleName, samFile, sort_order = "coordinate", validation_stringency = "LENIENT", create_index = True, clobber = False, outputDirectory = ""):
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
        self.validation_stringency = validation_stringency
        self.create_index = create_index
        self.sort_order = sort_order
        self.clobber = clobber
        #SANITY TEST ALL THE THINGS
        checkForRequiredFile(self.samFile, "SAM file to convert to BAM")
        if not type(self.create_index) == bool:
            raise RuntimeError("Create index value should be passed as a boolean. Passed value: %s" %(self.create_index))
        #DONE SANITY CHECKING. FOR NOW.
        if self.create_index:
            self.create_index = "true"
        else:
            self.create_index = "false"
        self.makeAndCheckOutputFileNames()
        self.samToBamCommand = self.createPicardCommand()
        
    def makeAndCheckOutputFileNames(self):
        self.bamOut = self.outputDirectory + self.sampleName + ".bam"
        self.clobber = checkForOverwriteRisk(self.bamOut, self.sampleName, self.clobber)
        
    def createPicardCommand(self):
        import genericRunners
        flagValues = {"I" : self.samFile,
                      "O" : self.bamOut,
                      "VALIDATION_STRINGENCY" : self.validation_stringency,
                      "CREATE_INDEX" : self.create_index,
                      "SORT_ORDER" : self.sort_order}
        picardArgs = [programPaths["java"], "-Xmx1g", "-jar", programPaths["picard"], "SortSam", flagValues]
        argumentFormatter = genericRunners.ArgumentFormatter(picardArgs, "=")
        picardCommand = argumentFormatter.argumentString
        return picardCommand
    
class Deduplicate(object):
    
    def __init__(self, sampleName, bamIn, validation_stringency = "LENIENT", create_index = True, clobber = False, outputDirectory = ""):
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
        self.bamIn = bamIn
        self.validation_stringency = validation_stringency
        self.create_index = create_index
        self.clobber = clobber
        #SANITY TEST ALL THE THINGS
        checkForRequiredFile(self.bamIn, "BAM file to deduplicate")
        if not type(self.create_index) == bool:
            raise RuntimeError("Create index value should be passed as a boolean. Passed value: %s" %(self.create_index))
        #DONE SANITY CHECKING. FOR NOW.
        if self.create_index:
            self.create_index = "true"
        else:
            self.create_index = "false"
        self.makeAndCheckOutputFileNames()
        self.deduplicateCommand = self.createPicardCommand()
        
    def makeAndCheckOutputFileNames(self):
        self.bamOut = self.outputDirectory + stripDirectoryAndExtension(self.bamIn) + ".deduped.bam"
        self.clobber = checkForOverwriteRisk(self.bamOut, self.sampleName, self.clobber)
        self.metricsOut = self.outputDirectory + stripDirectoryAndExtension(self.bamIn) + "dedupe.metrics"
        self.clobber = checkForOverwriteRisk(self.metricsOut, self.sampleName, self.clobber)
        
    def createPicardCommand(self):
        import genericRunners
        flagValues = {"I" : self.bamIn,
                      "O" : self.bamOut,
                      "VALIDATION_STRINGENCY" : self.validation_stringency,
                      "CREATE_INDEX" : self.create_index,
                      "METRICS_FILE" : self.metricsOut}
        picardArgs = [programPaths["java"], "-Xmx1g", "-jar", programPaths["picard"], "MarkDuplicates", flagValues]
        argumentFormatter = genericRunners.ArgumentFormatter(picardArgs, "=")
        picardCommand = argumentFormatter.argumentString
        return picardCommand
    
class MPileup(object):
    
    def __init__(self, sampleName, bamFile, refGenomeFasta, disablePerBaseAlignmentQuality = True, minBaseQuality = 20, maxDepth = 1000000, countOrphans = True, clobber = False, outputDirectory = ""):
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
        checkForRequiredFile(self.bamFile, "BAM file to analyze")
        if not type(self.refGenomeFasta) == str:
            raise RuntimeError("Reference genome fasta file name should be a string. Passed %s" %(self.refGenomeFasta))
        checkForRequiredFile(self.refGenomeFasta, "reference genome fasta")
        checkForRequiredFile(self.refGenomeFasta + ".fai", "reference genome fasta index file", "Please move the index file to this location or create one using samtools faidx.")
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
        self.mPileupOut = self.outputDirectory + stripDirectoryAndExtension(self.bamFile) + ".mpileup"
        self.clobber = checkForOverwriteRisk(self.mPileupOut, self.sampleName, self.clobber)
        
    def createSamtoolsCommand(self):
        import genericRunners
        flagValues = {"-B" : self.disablePerBaseAlignmentQuality,
                      "-Q" : self.minBaseQuality,
                      "-d" : self.maxDepth,
                      "-A" : self.countOrphans,
                      "-f" : self.refGenomeFasta}
        samtoolsArgs = [programPaths["samtools"], "mpileup", flagValues, self.bamFile, ">", self.mPileupOut]
        argumentFormatter = genericRunners.ArgumentFormatter(samtoolsArgs)
        samtoolsCommand = argumentFormatter.argumentString
        return samtoolsCommand
    
class ExtractVariantsTumor(object):
    
    def __init__(self, sampleName, pileupInput, clobber = False, minSupport = 0, requireDoubleStranded = False, outputDirectory = ""):
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
        checkForRequiredFile(self.pileupInput, "VCF input file")
        #done sanity checking
        self.makeAndCheckOutputFileNames()
        self.extractCommand = self.createPileupToVCFCommand()    
    
    def makeAndCheckOutputFileNames(self):
        self.variantsOut = self.outputDirectory + stripDirectoryAndExtension(self.pileupInput) + ".variants"
        self.targetList = self.outputDirectory + stripDirectoryAndExtension(self.pileupInput) + ".targets"
        self.clobber = checkForOverwriteRisk(self.variantsOut, self.sampleName, self.clobber)
        self.clobber = checkForOverwriteRisk(self.targetList, self.sampleName, self.clobber)
    
    def createPileupToVCFCommand(self):
        import genericRunners
        flagValues = {"-f" : self.pileupInput,
                      "-o" : self.variantsOut,
                      "-n" : self.minSupport,
                      "-t" : self.targetList,
                      "-d" : self.requireDoubleStranded}
        pileupCommandArgs = [programPaths["python3"], programPaths["extractVariants"], flagValues]
        argumentFormatter = genericRunners.ArgumentFormatter(pileupCommandArgs)
        pileupCommand = argumentFormatter.argumentString
        return pileupCommand
    
class ExtractVariantsNormal(object):
    
    def __init__(self, sampleName, pileupInput, targetList, comparison, minSupport = 0, clobber = False, outputDirectory = ""):
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
        checkForRequiredFile(self.pileupInput, "VCF input file")
        checkForRequiredFile(self.targetList, "List of targets from tumor to do ")
        #done sanity checking
        self.makeAndCheckOutputFileNames()
        self.extractCommand = self.createPileupToVCFCommand()
    
    def makeAndCheckOutputFileNames(self):
        self.variantsOut = self.outputDirectory + stripDirectoryAndExtension(self.pileupInput) + self.comparison + ".variants"
        self.clobber = checkForOverwriteRisk(self.variantsOut, self.sampleName, self.clobber)
    
    def createPileupToVCFCommand(self):
        import genericRunners
        flagValues = {"-f" : self.pileupInput,
                      "-o" : self.variantsOut,
                      "-n" : self.minSupport,
                      "-m" : self.targetList}
        pileupCommandArgs = [programPaths["python3"], programPaths["extractVariants"], flagValues]
        argumentFormatter = genericRunners.ArgumentFormatter(pileupCommandArgs)
        pileupCommand = argumentFormatter.argumentString
        return pileupCommand
    
class CombineExtractedVariants(object):
    
    def __init__(self, sampleName, tumorFileName, normalFileName, comparison, clobber = False, outputDirectory = ""):
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
        checkForRequiredFile(self.tumorFile, "Tumor input file")
        checkForRequiredFile(self.normalFile, "Normal input file")
        #done sanity checking
        self.makeAndCheckOutputFileNames()
        self.combineCommand = self.createPileupToVCFCommand()    
    
    def makeAndCheckOutputFileNames(self):
        self.vcfOut = self.outputDirectory + stripDirectoryAndExtension(self.normalFile) + self.comparison + ".vcf"
        self.clobber = checkForOverwriteRisk(self.vcfOut, self.sampleName, self.clobber)
    
    def createPileupToVCFCommand(self):
        import genericRunners
        flagValues = {"-t" : self.tumorFile,
                      "-n" : self.normalFile,
                      "-o" : self.vcfOut}
        combineCommandArgs = [programPaths["python3"], programPaths["combineVariants"], flagValues]
        argumentFormatter = genericRunners.ArgumentFormatter(combineCommandArgs)
        combineCommand = argumentFormatter.argumentString
        return combineCommand
        
        