#!/usr/bin/env python3

import os

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
    fileName = stripDirectory(fileName)
    fileName = stripExtension(fileName)
    return fileName

def baseFileName(fileName):
    import os
    return fileName.split(os.sep)[-1].split(".")[0]

def directory(fileName):
    import os
    return os.sep.join(fileName.split(os.sep)[:-1])

def stripExtension(fileName, zipDoubleExtension = True):
    zipExtensions = ["gz", "bz", "zip"]
    fileName = fileName.split(".")
    if len(fileName) < 3 or not fileName[-1] in zipExtensions or not zipDoubleExtension:
        return ".".join(fileName[:-1])
    else:
        return ".".join(fileName[:-2])

def stripDirectory(fileName):
    import os
    return fileName.split(os.sep)[-1]

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
    
def checkForRequiredDirectory(directoryName, directoryDescription, instructions = ""):
    import os
    if not os.path.isdir(directoryDescription):
        if instructions:
            if not instructions.startswith(" "):
                instructions = " " + instructions
        raise FileNotFoundError("Unable to find %s.%s Expected file: %s" %(directoryDescription, instructions, directoryName))
    
def checkTypes(values, requiredTypes):
    if not type(values) == list:
        values = [values]
    typesList = [type(value) for value in values]
    if not type(requiredTypes) in (list, tuple):
        requiredTypes = [requiredTypes]
    for value in values:
        assert type(value) in requiredTypes, "Incorrect value type found while type checking. Required: %s\nFound:\n%s\%s" %(requiredTypes, values, typesList)

class ArgumentFormatter(object):
    
    def __init__(self, arguments, keyValues = False, delimiter = ","):
        self.arguments = arguments
        self.delimiter = delimiter
        self.argumentList = self.formatArguments(tuple(self.arguments), keyValues, self.delimiter)
        self.argumentString = self.formatArgumentString(self.argumentList)
    
    def formatArguments(self, arguments, keyValues = False, delimiter = ","):  #doing the tuple and untuple thing to lower chance of having the original changed
        if keyValues:
            if type(keyValues) == bool:
                raise RuntimeError("If using key value outputs, a joining character (such as an equal sign or colon) must be passed.")
        arguments = list(arguments)
        argumentBlockList = []
        for argument in arguments:
            if type(argument) == str:
                argumentBlockList += [argument]
            elif type(argument) == list:
                argumentStringList = [str(item) for item in argument]
                argumentBlockList.append(delimiter.join(argumentStringList))
            elif type(argument) == dict:
                dictList = []
                keys = list(argument.keys())
                for key in keys:
                    directAdd = False
                    dictList += [key]
                    if type(argument[key]) == str:
                        addString = argument[key]
                    elif type(argument[key]) == list:
                        if key.upper().startswith("FLAGGEDLIST"):
                            directAdd = True
                            del dictList[-1]
                            flaggedList = []
                            if keyValues or len(argument[key]) == 3:
                                spacer = keyValues  #if keyValues was not set, or even if we need to override it, we will do so immediately below
                                if len(argument[key]) == 3:  #allows the user to pass an additional item in the list to override any spacer
                                    spacer = argument[key][2]
                            else:
                                spacer = " "
                            for item in argument[key][1]:
                                flaggedList.append(argument[key][0] + spacer + item)
                            addString = " ".join(flaggedList)
                        else:
                            argument[key] = [str(item) for item in argument[key]]
                            addString = delimiter.join(argument[key])
                    elif type(argument[key]) in [float, int]:
                        addString = str(argument[key])
                    elif type(argument[key]) == bool:
                        addString = ""
                        if not argument[key]:
                            del dictList[-1]  #if it was false on a boolean argument, remove the flag we put in at the start of this method.
                    if keyValues and not directAdd:
                        dictList[-1] += keyValues + addString
                    else:
                        if addString:
                            dictList.append(addString)
                argumentBlockList += dictList
        return argumentBlockList                
            
    def formatArgumentString(self, arguments):
        return " ".join(arguments)
    
    def __str__(self):
        return self.argumentString
    
class WorkflowReturn(object):
    
    def __init__(self, data, jobID, clobber, intermediateData = {}, intermediateJobs = {}):
        self.data = data
        self.jobID = jobID
        self.clobber = clobber
        self.intermediateData = intermediateData
        self.intermediateJobs = intermediateJobs
        checkTypes(self.data, [str, list])
        checkTypes(self.clobber, [bool, type(None)])
        checkTypes(self.jobID, [int, list])
        
    def addJob(self, workflowReturnIn):
        checkTypes(workflowReturnIn.data, [str, list])
        checkTypes(workflowReturnIn.jobID, [list, int])
        checkTypes(workflowReturnIn.clobber, [bool, type(None)])
        if type(self.jobID) == int:
            self.jobID = [self.jobID]
        if type(self.data) == str:
            self.data = [self.data]
        if type(workflowReturnIn.data) == str:
            workflowReturnIn.data = [workflowReturnIn.data]
        self.data += (workflowReturnIn.data)
        if type(workflowReturnIn.jobID) == int:
            workflowReturnIn.jobID = [workflowReturnIn.jobID]
        self.jobID += (workflowReturnIn.jobID)
        if self.clobber == None:  #if no clobber value has been firmly set, override it with whatever comes in, even if it's staying as none
            self.clobber = workflowReturnIn.clobber
        elif not self.clobber and workflowReturnIn.clobber: #if self.clobber is false and the incoming value is true, override it with the new value (this probably won't happen)
            self.clobber = clobber
            
    def add(self, data, jobID, clobber):
        checkTypes(data, [str, list])
        checkTypes(jobID, [int, list])
        checkTypes(clobber, [bool, type(None)])
        if type(self.jobID) == int:
            self.jobID = [self.jobID]
        if type(jobID) == int:
            jobID = [jobID]
        if type(self.data) == str:
            self.data = [self.data]
        if type(data) == str:
            data = [data]
        self.data += (data)
        self.jobID += (jobID)
        if self.clobber == None:  #if no clobber value has been firmly set, override it with whatever comes in, even if it's staying as none
            self.clobber = clobber
        elif not self.clobber and clobber: #if self.clobber is false and the incoming value is true, override it with the new value (this probably won't happen)
            self.clobber = clobber
        
class EmptyWorkflowReturn(WorkflowReturn):
    
    def __init__(self):
        self.data = ""
        self.jobID = []
        self.clobber = None
        self.intermediateData = {}
        self.intermediateJobs = {}
        checkTypes(self.data, [str, list])
        checkTypes(self.clobber, [bool, type(None)])
        checkTypes(self.jobID, [int, list])
        
    def __bool__(self):
        return False
        
class JobList(object):
    
    def __init__(self, tempDir):
        import os
        import pickle
        self.tempDir = tempDir
        if not self.tempDir.endswith(os.sep):
            self.tempDir = self.tempDir + os.sep
        self.fileName = self.tempDir + "jobs.pkl"
        self.file = open(self.fileName, 'w+b') #need to read and write as a binary
        placeHolder = ["Zeroth job placeholder"] #need a placeholder for the zeroth position because the scheduler is indexed to job 1 and can't handle a job zero
        pickle.dump(placeHolder, self.file)
        self.file.seek(0)  #return to the start of the file for when we need to rewrite it
        
    def addJob(self, command):
        import pickle
        self.file.seek(0)  #reset to the start of the file
        jobList = pickle.load(self.file)
        self.file.seek(0)  #reset for an impending overwrite
        newCommandIndex = len(jobList) #this will be the last index + 1, which will be the index of the next job
        jobList.append(command)
        pickle.dump(jobList, self.file)
        self.file.seek(0)  #reset file for whatever we do next
        return newCommandIndex  #return the index number of the new job
    
    def close(self):
        self.file.close()
        
def createTempDir(subjectID, useDir = False):
    import re
    import os
    import datetime
    successful = False
    if not useDir:
        useDir = os.getcwd()
    while not successful:
        currenttime = datetime.datetime.now()
        currenttime = str(currenttime)
        currenttime = re.sub(r'\W','',currenttime)
        tempDir = useDir + os.sep + "." + subjectID + currenttime
        if os.path.isdir(tempDir):
            continue
        try:
            os.mkdir(tempDir)
        except OSError:
            continue
        successful = True
    os.mkdir(tempDir + "/outputs")
    return tempDir
    
def calculateCoresFromFileSize(files, gigsPerCore, maxCores = 8):
    import os
    if type(files) == str:
        files = [files]
    totalSize = 0
    for file in files:
        if file:
            totalSize += os.path.getsize(file)
            #print(file + " : " + str(os.path.getsize(file)))
    cores = int(-(-totalSize // (gigsPerCore * 1000000000)))
    cores = min([cores, maxCores])
    #print("Cores: %s" %(cores))
    if cores < 1:
        return 1
    else:
        return cores
    