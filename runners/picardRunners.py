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
                "bam-readcount" : os.getcwd() + "bin/bam-readcount/bin/bam-readcount"}
global picardPath
picardPath = os.getcwd() + "/bin/picard-tools-2.1.1/picard.jar"

class sortSAMtoBAM(object):
    
    def __init__(self, sampleName, samFile, sort_order = "coordinate", validation_stringency = "LENIENT", create_index = True, clobber = False, outputDirectory = ""):
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
        self.samFile = samFile
        self.validation_stringency = validation_stringency
        self.create_index = create_index
        self.sort_order = sort_order
        self.clobber = clobber
        #SANITY TEST ALL THE THINGS
        runnerSupport.checkForRequiredFile(self.samFile, "SAM file to convert to BAM")
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
        import runnerSupport
        self.bamOut = self.outputDirectory + self.sampleName + ".bam"
        self.clobber = runnerSupport.checkForOverwriteRisk(self.bamOut, self.sampleName, self.clobber)
        
    def createPicardCommand(self):
        import runnerSupport
        flagValues = {"I" : self.samFile,
                      "O" : self.bamOut,
                      "VALIDATION_STRINGENCY" : self.validation_stringency,
                      "CREATE_INDEX" : self.create_index,
                      "SORT_ORDER" : self.sort_order}
        picardArgs = [programPaths["java"], "-Xmx1g", "-jar", picardPath, "SortSam", flagValues]
        argumentFormatter = runnerSupport.ArgumentFormatter(picardArgs, "=")
        picardCommand = argumentFormatter.argumentString
        return picardCommand
    
class MergeSAMFiles(object):
    
    def __init__(self, sampleName, bamFiles, sort_order = "coordinate", validation_stringency = "LENIENT", create_index = True, clobber = False, outputDirectory = ""):
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
        self.bamFiles = bamFiles
        self.validation_stringency = validation_stringency
        self.create_index = create_index
        self.sort_order = sort_order
        self.clobber = clobber
        #SANITY TEST ALL THE THINGS
        for file in self.bamFiles:
            runnerSupport.checkForRequiredFile(file, "BAM file to merge")
        if not type(self.create_index) == bool:
            raise RuntimeError("Create index value should be passed as a boolean. Passed value: %s" %(self.create_index))
        #DONE SANITY CHECKING. FOR NOW.
        if self.create_index:
            self.create_index = "true"
        else:
            self.create_index = "false"
        self.makeAndCheckOutputFileNames()
        self.mergeSAMCommand = self.createPicardCommand()
        
    def makeAndCheckOutputFileNames(self):
        import runnerSupport
        self.bamOut = self.outputDirectory + self.sampleName + ".merged.bam"
        self.clobber = runnerSupport.checkForOverwriteRisk(self.bamOut, self.sampleName, self.clobber)
        
    def createPicardCommand(self):
        import runnerSupport
        flagValues = {"flaggedlist" : ["I", self.bamFiles],
                      "O" : self.bamOut,
                      "VALIDATION_STRINGENCY" : self.validation_stringency,
                      "CREATE_INDEX" : self.create_index,
                      "SORT_ORDER" : self.sort_order}
        picardArgs = [programPaths["java"], "-Xmx1g", "-jar", picardPath, "MergeSamFiles", flagValues]
        argumentFormatter = runnerSupport.ArgumentFormatter(picardArgs, "=")
        picardCommand = argumentFormatter.argumentString
        return picardCommand
    
class AddReadGroups(object):
    
    def __init__(self, sampleName, bamIn, sort_order = "coordinate", rgid = "1", rglb = "defaultLibrary", rgpl = "Illumina", rgpu = "defaultBarcode", rgsm = False, validation_stringency = "LENIENT", create_index = True, clobber = False, outputDirectory = ""):
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
        self.bamIn = bamIn
        self.validation_stringency = validation_stringency
        self.create_index = create_index
        self.sort_order = sort_order
        self.clobber = clobber
        #SANITY TEST ALL THE THINGS
        runnerSupport.checkForRequiredFile(self.bamIn, "BAM file for adding readgroups")
        if not type(self.create_index) == bool:
            raise RuntimeError("Create index value should be passed as a boolean. Passed value: %s" %(self.create_index))
        #DONE SANITY CHECKING. FOR NOW.
        # if self.create_index:
        #     self.create_index = "true"
        # else:
        #     self.create_index = "false"
        if not rgid:
            raise RuntimeError("RGID value must be given.")
        self.rgid = str(rgid)
        if not rglb:
            raise RuntimeError("RGLB value must be given.")
        self.rglb = str(rglb)
        if not rgpl:
            raise RuntimeError("RGLB value must be given.")
        self.rgpl = str(rgpl)
        if not rgpu:
            raise RuntimeError("RGPU value must be given.")
        self.rgpu = str(rgpu)
        if not rgsm:
            self.rgsm = str(self.sampleName)
        else:
            self.rgsm = str(rgsm)
        if self.create_index:
            self.create_index = "true"
        else:
            self.create_index = "false"
        self.makeAndCheckOutputFileNames()
        self.addReadGroupsCommand = self.createPicardCommand()
        
    def makeAndCheckOutputFileNames(self):
        import runnerSupport
        self.bamOut = self.outputDirectory + runnerSupport.stripDirectoryAndExtension(self.bamIn) + ".rdgp" + ".bam"
        self.clobber = runnerSupport.checkForOverwriteRisk(self.bamOut, self.sampleName, self.clobber)
        
    def createPicardCommand(self):
        import runnerSupport
        flagValues = {"I" : self.bamIn,
                      "O" : self.bamOut,
                      "VALIDATION_STRINGENCY" : self.validation_stringency,
                      "CREATE_INDEX" : self.create_index,
                      "SORT_ORDER" : self.sort_order,
                      "RGID" : self.rgid,
                      "RGPL" : self.rgpl,
                      "RGPU" : self.rgpu,
                      "RGSM" : self.rgsm,
                      "RGLB" : self.rglb}
        picardArgs = [programPaths["java"], "-Xmx1g", "-jar", picardPath, "AddOrReplaceReadGroups", flagValues]
        argumentFormatter = runnerSupport.ArgumentFormatter(picardArgs, "=")
        picardCommand = argumentFormatter.argumentString
        return picardCommand
    
class Deduplicate(object):
    
    def __init__(self, sampleName, bamIn, validation_stringency = "LENIENT", create_index = True, revertToBaseName = False, clobber = False, outputDirectory = ""):
        import runnerSupport
        self.revertToBaseName = revertToBaseName
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
        runnerSupport.checkForRequiredFile(self.bamIn, "BAM file to deduplicate")
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
        import runnerSupport
        if self.revertToBaseName:
            self.bamOut = self.outputDirectory + self.sampleName + ".deduped.bam"
            self.metricsOut = self.outputDirectory + self.sampleName + ".dedupe.metrics"
        else:
            self.bamOut = self.outputDirectory + runnerSupport.stripDirectoryAndExtension(self.bamIn) + ".deduped.bam"
            self.metricsOut = self.outputDirectory + runnerSupport.stripDirectoryAndExtension(self.bamIn) + ".dedupe.metrics"
        self.clobber = runnerSupport.checkForOverwriteRisk(self.bamOut, self.sampleName, self.clobber)
        self.clobber = runnerSupport.checkForOverwriteRisk(self.metricsOut, self.sampleName, self.clobber)
        
    def createPicardCommand(self):
        import runnerSupport
        flagValues = {"I" : self.bamIn,
                      "O" : self.bamOut,
                      "VALIDATION_STRINGENCY" : self.validation_stringency,
                      "CREATE_INDEX" : self.create_index,
                      "METRICS_FILE" : self.metricsOut}
        picardArgs = [programPaths["java"], "-Xmx1g", "-jar", picardPath, "MarkDuplicates", flagValues]
        argumentFormatter = runnerSupport.ArgumentFormatter(picardArgs, "=")
        picardCommand = argumentFormatter.argumentString
        return picardCommand
    
class BuildBAMIndex(object):
    
    def __init__(self, bamFile, validation_stringency = "LENIENT", create_index = True):
        import runnerSupport
        self.bamFile = bamFile
        self.validation_stringency = validation_stringency
        #SANITY TEST ALL THE THINGS
        runnerSupport.checkForRequiredFile(self.bamFile, "BAM file to index")
        #DONE SANITY CHECKING. FOR NOW.
        self.indexCommand = self.createPicardCommand()
        #not checking for future file collisions here.  if there is a collision, we want it to be current.
        
    def createPicardCommand(self):
        import runnerSupport
        flagValues = {"I" : self.bamFile,
                      "VALIDATION_STRINGENCY" : self.validation_stringency}
        picardArgs = [programPaths["java"], "-Xmx1g", "-jar", picardPath, "BuildBamIndex", flagValues]
        argumentFormatter = runnerSupport.ArgumentFormatter(picardArgs, "=")
        picardCommand = argumentFormatter.argumentString
        return picardCommand
    