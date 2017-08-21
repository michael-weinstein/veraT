#!/usr/bin/env python3

import os
global runnerRoot
runnerRoot = os.sep.join(os.path.abspath(__file__).split(os.sep)[:-3]) + os.sep
global programPaths
programPaths = {"netMHC" : runnerRoot + "/bin/netMHC-4.0/netMHC",
                "netMHCHome" : runnerRoot + "/bin/netMHC-4.0/",
                "python3" : "/u/local/apps/python/3.4.3/bin/python3",
                "netMHCRunner" : os.path.abspath(__file__),
                "permittedHLAList" : runnerRoot + "/bin/netMHC-4.0/data/allelelist"}

class CheckArgs():  #class that checks arguments and ultimately returns a validated set of arguments to the main program
    
    def __init__(self):
        import argparse
        import os
        import re
        parser = argparse.ArgumentParser()
        parser.add_argument("-i", "--inputVariantFile", help = "Pickle of variants with added HLA data", required = True)
        parser.add_argument("-o", "--output", help = "Output file base name", required = True)
        parser.add_argument("-t", "--tempdir", help = "Temporary directory for netMHC", required = True)
        rawArgs = parser.parse_args()
        inputVariantFile = rawArgs.inputVariantFile
        if not os.path.isfile(inputVariantFile):
            raise FileNotFoundError("Unable to find input pickle file at %s" %inputVariantFile)
        self.inputVariantFile = inputVariantFile
        self.output = rawArgs.output
        tempdir = rawArgs.tempdir
        if not os.path.exists(tempdir):
            try:
                os.mkdir(tempdir)
            except Exception as e:
                raise RuntimeError("Error making temp directory for netMHC. Attempt returned %s" %e)
        else:
            if not os.path.isdir(tempdir):
                raise RuntimeError("File already exists at given temporary directory location. Unable to create %s" %tempdir)
        self.tempdir = tempdir

def getAlleleList():
    alleleFile = open(programPaths["permittedHLAList"], 'r')
    alleleList = []
    for line in alleleFile:
        line = line.strip()
        if not line:
            continue
        alleleList.append(line.split()[0])
    return alleleList
        
def getFormattedHLAList(hlaDict):
    def processHLA(unformattedHLA):
        import re
        formattedHLA = unformattedHLA
        formattedHLA = formattedHLA.upper()
        formattedHLA = formattedHLA.replace("HLA","")
        formattedHLA = re.sub("\W", "", formattedHLA)
        if not re.match("^[ABC]\d\d\d\d\d?$", formattedHLA):
            raise RuntimeError("Unable to read HLA given: %s" %unformattedHLA)
        return "HLA-" + formattedHLA
    hlaMolecules = list(hlaDict.keys())
    formattedHLAList = []
    for molecule in hlaMolecules:
        if hlaDict[molecule]:
            for genotype in hlaDict[molecule]:
                formattedHLAList.append(processHLA(genotype))
    return formattedHLAList
        
def setEnvironmentVariables(tempdir):
    import os
    os.environ['NMHOME'] = programPaths["netMHCHome"]
    os.environ['TMPDIR'] = tempdir
    
def makeTempDir(containingDirectory):
    import os
    import datetime
    import re
    created = False
    attempts = 0
    while not created and attempts < 10:
        tempdirFolder = "netMHCworking%s" %(datetime.datetime.now())
        tempdirFolder = re.sub("\W", "", tempdirFolder)
        tempdir = containingDirectory + os.sep + tempdirFolder
        try:
            os.mkdir(tempdir)
            created = True
        except:
            attempts += 1
    if not created:
        tempdir = containingDirectory + os.sep + "netMHCworking%s" %(datetime.datetime.now())
        os.mkdir(tempdir)  #if it works, fine, but if it's failing, it will throw the exception here
    return tempdir
    
class netMHCWrapper(object):
    
    def __init__(self, sampleName, variantFile, clobber = False, outputDirectory = ""):
        import runnerSupport
        import os
        self.variantFile = variantFile
        self.sampleName = sampleName
        self.clobber = clobber
        self.peptideFileArgList = []
        peptideFileListForValidation = []
        self.netMHCtempdir = makeTempDir(outputDirectory)
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
        runnerSupport.checkForRequiredFile(self.variantFile, "Input variant file with HLA data")
        for file in peptideFileListForValidation:
            runnerSupport.checkForRequiredFile(file, "Peptide list file")
        assert os.path.isdir(self.netMHCtempdir), "Unable to find netMHC working folder directory %s" %self.netMHCtempdir
        #DONE SANITY CHECKING. FOR NOW.
        self.makeAndCheckOutputFileNames()
        self.netMHCWrapperCommand = self.createNetMHCWrapperCommand()
    
    def makeAndCheckOutputFileNames(self):
        import runnerSupport
        self.rawNetMHCOutput = self.outputDirectory + self.sampleName + ".rawNetMHC.predictions"
        self.peptideTableOut = self.outputDirectory + self.sampleName + ".peptidePredictions.txt"
        self.clobber = runnerSupport.checkForOverwriteRisk(self.rawNetMHCOutput, self.sampleName, self.clobber)
        self.baseSampleName = self.outputDirectory + self.sampleName

    def createNetMHCWrapperCommand(self):
        import runnerSupport
        flagValues = {"-i" : self.variantFile,
                      "-o" : self.baseSampleName,
                      "-t" : self.netMHCtempdir}
        args = [programPaths["python3"], programPaths["netMHCRunner"], flagValues]
        argumentFormatter = runnerSupport.ArgumentFormatter(args)
        netMHCWrapperCommand = argumentFormatter.argumentString
        return netMHCWrapperCommand

class PeptideListMakerArgs(object):
    
    def __init__(self, variantsFile, output):
        self.variantsFile = variantsFile
        self.output = output

def loadVariantPickle(variantPickleFile):
    import pickle
    file = open(variantPickleFile, 'rb')
    data = pickle.load(file)
    file.close()
    return data

def formatDictionaryForOutput(dictionary, delimiter = ", "):
    keys = sorted(list(dictionary.keys()))
    outputList = []
    for key in keys:
        outputList.append(key + ":" + str(dictionary[key]))
    return delimiter.join(outputList)

def createOutputTable(variantPickle):
    import numpy
    import variantDataHandler
    commentLines = []
    hlaCalls = []
    for molecule in sorted(list(variantPickle["hla"].keys())):
        hlaCalls += sorted(variantPickle["hla"][molecule])
    hlaLine = "HLA Calls:" + ", ".join(hlaCalls)
    commentLines.append(hlaLine)
    commentLines.append("Mutation class counts: " + formatDictionaryForOutput(variantPickle["fused"]["oncotatorClassCounts"]))
    commentLines.append("Mutation type counts: " + formatDictionaryForOutput(variantPickle["fused"]["oncotatorTypeCounts"]))
    usingRNA = "RNASupport" in variantPickle
    headerLine = ["#chrom", "pos", "ref", "alt", "variantType", "variantClass", "gene", "description",  "genomeChange", "transcriptID", "cDNAChange", "exon", "proteinChange", "HLA", "altPeptide",  "alt1-log50k(aff)", "altAffinity(nM)", "altPercentRank", "refPeptide", "ref1-log50k(aff)", "refAffinity", "refPercentRank", "logRatio", "affinityRatio", "deltaPercent", "tumorSupportPercent", "tumorSupportReads", "tumorDepth", "normalSupportPercent", "normalSupportReads", "normalDepth"]
    rnaColumns = ["tumorRNASupportScore", "tumorRNASupportPercent", "tumorRNASupport", "tumorRNADepth", "tumorRNApValue", "tumorRNAOddsRatio"]
    if usingRNA:
        headerLine += rnaColumns
    locusList = list(variantPickle["fused"]["variants"].keys())
    variantDataHandler.sortVariantDataTuples(locusList)
    outputData = []
    for locus in locusList:
        chromosome, position, ref, alt = locus
        oncotatorLocus = (chromosome.replace("chr",""), position, ref, alt)
        if not oncotatorLocus in variantPickle["fused"]["oncotatorData"]:
            #print("Could not find %s in oncotator data" %(str(oncotatorLocus)))
            continue
        #else:
            #print("Found %s in oncotatorData" %(str(oncotatorLocus)))
        oncotatorData = variantPickle["fused"]["oncotatorData"][oncotatorLocus]
        dnaReadData = variantPickle["fused"]["variants"][locus]
        mutantPeptideList = variantPickle["fused"]["neoepitopes"][oncotatorLocus]
        wildtypeProteinList = variantPickle["fused"]["wildtype"][oncotatorLocus]
        if usingRNA:
            rnaReadData = variantPickle["fused"]["RNASupport"][locus]
            try:
                tumorRNASupportPercent = rnaReadData.supportingReads/rnaReadData.totalDepth
            except ZeroDivisionError:
                tumorRNASupportPercent = numpy.nan
        for index, mutantPeptide in enumerate(mutantPeptideList):
            wildtypePeptide = wildtypeProteinList[index]
            mutantPredictionDict = variantPickle["fused"]["netMHC"][mutantPeptide]
            wildtypePredictionDict = variantPickle["fused"]["netMHC"][wildtypePeptide]
            hlaAlleles = sorted(list(mutantPredictionDict.keys()))
            for hla in hlaAlleles:
                mutantPrediction = mutantPredictionDict[hla]
                wildtypePrediction = wildtypePredictionDict[hla]
                try:
                    logRatio = mutantPrediction.log/wildtypePrediction.log
                except ZeroDivisionError:
                    logRatio = numpy.nan
                try:
                    affinityRatio = mutantPrediction.affinity/wildtypePrediction.affinity
                except ZeroDivisionError:
                    affinityRatio = numpy.nan
                deltaPercent = mutantPrediction.rank-wildtypePrediction.rank
                outputList = [chromosome, position, ref, alt, oncotatorData.variantType, oncotatorData.variantClassification, oncotatorData.gene, oncotatorData.description, oncotatorData.genomeChange, oncotatorData.transcript, oncotatorData.cDNAChange, oncotatorData.exon, oncotatorData.proteinChange, hla, mutantPeptide, mutantPrediction.log, mutantPrediction.affinity, mutantPrediction.rank, wildtypePeptide, wildtypePrediction.log, wildtypePrediction.affinity, wildtypePrediction.rank, logRatio, affinityRatio, deltaPercent, dnaReadData.tumorSupportingPercent, dnaReadData.tumorSupporting, dnaReadData.tumorDepth, dnaReadData.normalSupportingPercent, dnaReadData.normalSupporting, dnaReadData.normalDepth]
                if usingRNA:
                    rnaOutput = [rnaReadData.score, tumorRNASupportPercent, rnaReadData.supportingReads, rnaReadData.totalDepth, rnaReadData.pvalue, rnaReadData.oddsRatio]
                    outputList += rnaOutput
                outputData.append(outputList)
    return (commentLines, headerLine, outputData)

def formatOutputTableString(outputData, delimiter = "\t", newline = "\n"):
    formattedTable = [[str(item) for item in line] for line in outputData]
    formattedTable = [delimiter.join(line) for line in formattedTable]
    return newline.join(formattedTable)
        
def runNetMHC():
    import os
    import variantDataHandler
    import peptideListMaker
    import netMHCReader
    import sys
    args = CheckArgs()
    peptideListMakerArgs = PeptideListMakerArgs(args.inputVariantFile, args.output)
    peptideListDict = peptideListMaker.makePeptideList(peptideListMakerArgs)
    setEnvironmentVariables(args.tempdir)
    variantPickle = loadVariantPickle(args.inputVariantFile)
    hlaList = getFormattedHLAList(variantPickle["hla"])
    peptideLengthList = list(peptideListDict.keys())
    permittedAlleles = getAlleleList()
    reportedBadAlleles = []
    rawNetMHCOutput = args.output + ".rawNetMHC.predictions"
    teeCharacter = " > " #this will start a new file for the first output, then be changed to append the next after
    for length in peptideLengthList:
        for hla in hlaList:
            if not hla in permittedAlleles:
                if not hla in reportedBadAlleles:
                    print("[WARNING] %s is not an analyzable allele in netMHC. Skipping it." %hla)
                    reportedBadAlleles.append(hla)
                continue
            commandList = [programPaths["netMHC"], "-p", "-l %s" %length, "-a %s" %hla, peptideListDict[length], teeCharacter, rawNetMHCOutput]
            teeCharacter = " >> "  #setting this to append subsequent files (should already be this if not the first run)
            command = " ".join(commandList)
            print("RUNNING: " + command)
            returnCode = os.system(command)
            print("COMPLETED. Exit status %s" %returnCode)
            if not returnCode == 0:
                print("NetMHC exited with non-zero status: %s" %returnCode)
                sys.exit(1)
    variantPickle["fused"]["netMHC"] = netMHCReader.getPredictionTable(rawNetMHCOutput)
    commentLines, headerLine, outputData = createOutputTable(variantPickle)
    outputTableString = formatOutputTableString(outputData)
    headerLineString = "\t".join(headerLine)
    if commentLines:
        commentLines = ("##" + line for line in commentLines)
        commentLines = "\n".join(commentLines)
    outputFile = open(args.output + ".peptidePredictions.txt", 'w')
    if commentLines:
        outputFile.write(commentLines + "\n")
    outputFile.write(headerLineString + "\n")
    outputFile.write(outputTableString)
    outputFile.close()    

if __name__ == '__main__':
    runNetMHC()