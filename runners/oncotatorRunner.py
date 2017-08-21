#!/usr/bin/env python3

import os
global runnerRoot
runnerRoot = os.sep.join(os.path.abspath(__file__).split(os.sep)[:-2]) + os.sep
global programPaths
programPaths = {"oncotator" : runnerRoot + "/bin/oncotator/oncotatorInstallation/bin/oncotator",
                "python2" : "/u/local/apps/python/2.7.2/bin/python",
                "python3" : "/u/local/apps/python/3.4.3/bin/python3",
                "oncotatorVitualEnv" : runnerRoot + "/bin/oncotator/oncotatorVirtualEnv/bin/activate",
                "oncotatorDefaultDB" : runnerRoot + "/references/oncotatorDB/oncotator_v1_ds_Jan262014/",
                "python2LibPath" : "/u/local/apps/python/2.7.2/lib",
                "oncotatorDefaultGenome" : "hg19",
                "oncotatorRunner" : os.path.abspath(__file__)}

class CheckArgs():  #class that checks arguments and ultimately returns a validated set of arguments to the main program
    
    def __init__(self):
        import argparse
        import os
        parser = argparse.ArgumentParser()
        parser.add_argument("-i", "--inputVariantFile", help = "Input pickle of fused, sorted variants.", required = True)
        parser.add_argument("-o", "--oncotatorOutput", help = "Output file name", required = True)
        parser.add_argument("-d", "--oncotatorDB", help = "Oncotator database directory", default = programPaths["oncotatorDefaultDB"])
        parser.add_argument("-g", "--genome", help = "Genome", default = programPaths["oncotatorDefaultGenome"])
        parser.add_argument("-r", "--outputDirectory", help = "Output directory (needed for the log file)", required = True)
        rawArgs = parser.parse_args()
        inputVariantFile = rawArgs.inputVariantFile
        if not os.path.isfile(inputVariantFile):
            raise FileNotFoundError("Unable to find input pickle file at %s" %inputVariantFile)
        self.inputVariantFile = inputVariantFile
        self.oncotatorOutput = rawArgs.oncotatorOutput
        oncotatorDB = rawArgs.oncotatorDB
        if not os.path.isdir(oncotatorDB):
            raise FileNotFoundError("Unable to find oncotator DB directory: %s" %oncotatorDB)
        self.oncotatorDB = oncotatorDB
        self.genome = rawArgs.genome
        self.outputDirectory = rawArgs.outputDirectory

def setEnvironmentVariables():
    import os
    oldLD = os.environ['LD_LIBRARY_PATH']
    tempLD = oldLD.split(":")
    tempLD.append(programPaths["python2LibPath"])
    os.environ['LD_LIBRARY_PATH'] = ":".join(tempLD)
    os.environ['PYTHONPATH'] = runnerRoot + "/bin/oncotator/oncotatorInstallation/lib/python2.7/site-packages"
    
class OncotatorWrapper(object):
    
    def __init__(self, sampleName, variantFile, outputFile = False, dataBase = False, genome = False, clobber = False, outputDirectory = ""):
        import runnerSupport
        import os
        self.variantFile = variantFile
        self.clobber = clobber
        if genome:
            self.genome = genome
        else:
            self.genome = programPaths["oncotatorDefaultGenome"]
        self.outputFile = outputFile
        self.sampleName = sampleName
        if dataBase:
            self.dataBase = dataBase
        else:
            self.dataBase = programPaths["oncotatorDefaultDB"]
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
        runnerSupport.checkForRequiredFile(self.variantFile, "Input variant file for oncotator")
        assert os.path.isdir(self.dataBase), "Unable to find oncotator database directory %s" %self.dataBase
        #DONE SANITY CHECKING. FOR NOW.
        self.makeAndCheckOutputFileNames()
        self.oncotatorWrapperCommand = self.createOncotatorWrapperCommand()
        
    def makeAndCheckOutputFileNames(self):
        import runnerSupport
        if self.outputFile:
            self.oncotatorOut = self.outputFile
        else:
            self.oncotatorOut = self.outputDirectory + self.sampleName + ".oncotator.txt"
        self.clobber = runnerSupport.checkForOverwriteRisk(self.oncotatorOut, self.sampleName, self.clobber)

    def createOncotatorWrapperCommand(self):
        import runnerSupport
        flagValues = {"-i" : self.variantFile,
                      "-o" : self.oncotatorOut,
                      "-g" : self.genome,
                      "-d" : self.dataBase,
                      "-r" : self.outputDirectory}
        args = [programPaths["python3"], programPaths["oncotatorRunner"], flagValues]
        argumentFormatter = runnerSupport.ArgumentFormatter(args)
        oncotatorWrapperCommand = argumentFormatter.argumentString
        return oncotatorWrapperCommand

def runOncotator(variantFile, outputFile, dataBase, genome, outputDirectory):
    import os
    setEnvironmentVariables()
    activateVirtualEnvCommand = "source %s" %programPaths["oncotatorVitualEnv"]
    changeDirectoryCommand = "cd %s" %outputDirectory
    moduleLoaderCommand = ". /u/local/Modules/default/init/modules.sh; module load python/2.7"
    oncotatorCommand = "%s -v --db-dir %s %s %s %s" %(programPaths["oncotator"], dataBase, variantFile, outputFile, genome)
    fullCommand = "%s; %s; %s; %s" %(moduleLoaderCommand, activateVirtualEnvCommand, changeDirectoryCommand, oncotatorCommand)
    exitStatus = os.system(fullCommand)
    return exitStatus

if __name__ == '__main__':
    import sys
    args = CheckArgs()
    exitStatus = runOncotator(args.inputVariantFile, args.oncotatorOutput, args.oncotatorDB, args.genome, args.outputDirectory)
    print("Oncotator exited with status %s" %exitStatus)
    sys.exit(exitStatus)