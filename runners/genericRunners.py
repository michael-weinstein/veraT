#!/usr/bin/env python3

global genericRunnerPaths
genericRunnerPaths = {"python3" : "/u/local/apps/python/3.4.3/bin/python3",  #hoffman
                      #"python3" : "/Library/Frameworks/Python.framework/Versions/3.4/bin/python3"  #local for debugging
                      "qsub" : "/u/systems/UGE8.0.1vm/bin/lx-amd64/qsub",
                      "arrayWrapper" : "./arrayWrapper.py"}

class HoffmanJob(object):
    
    def __init__(self, dependencies, jobNumber, jobName, tempDir, emailAddress, emailConditions, cores = 1, memory = 6, maxRetries = 10, mock = False):
        import os
        self.dependencies = dependencies
        if type(self.dependencies) in (int, str):
            self.dependencies = [self.dependencies]
        self.jobName = jobName
        self. jobNumber = jobNumber
        self.emailAddress = emailAddress
        self.emailConditions = emailConditions
        self.mock = mock
        if not type(self.jobNumber) == int:
            raise RuntimeError("Job number must be an integer")
        self.cores = cores
        if not type(cores) == int:
            raise RuntimeError("Core count must be an integer")
        self.memory = memory
        if not type(memory) == int or type(memory) == float:
            raise RuntimeError("Memory limit must be a number")
        self.memoryPerCore = self.calculateMemoryPerCore()
        self.tempDir = tempDir
        if not os.path.isdir(self.tempDir):
            raise RuntimeError("Unable to find working directory")
        self.jobID = self.submitJob(maxRetries)
        
    def calculateMemoryPerCore(self):  #because pe shared mode tends to give you your desired memory per core
        if self.cores == 1:
            return self.memory
        if self.memory % self.cores == 0:  #checking to see if it is evenly divisible
            return self.memory // self.cores
        return round(float(self.memory)/self.cores, 1)
        
    # def clearedDependencies(self, dependencies):
    #     import os
    #     for dependency in dependencies:
    #         if not os.path.isfile(dependency):
    #             return False
    #         if os.path.getsize(dependency) == 0:  #make sure it's not an empty file
    #             return False
    #     return True
    # 
    # def pulledFromJobBoard(self, jobNumber):  #job board will essentially have files as job tokens.  Only one process should be able to delete each file.
    #     import os
    #     jobBoardFile = self.tempDir + os.sep + "jobBoard" + os.sep + str(jobNumber) + ".job"
    #     try:
    #         os.remove(jobBoardFile)
    #     except FileNotFoundError:
    #         return False
    #     return True
    # 
    # def jobReady(self):
    #     if self.clearedDependencies(self.dependencies):
    #         if self.pulledFromJobBoard(self.jobNumber):  #very important they get tested in this order.  This line will take the token.
    #            return True
            
    def submitJob(self, maxRetries = 10):
        import os
        import subprocess
        peArg = []
        if self.cores > 1:
            peArg = ["-pe", "shared", str(self.cores)]
        limitsArgs = ["-l", "h_rt=23:59:59,h_data=" + str(self.memory) + "G"]
        jobRangeArg = ["-t", str(self.jobNumber) + "-" + str(self.jobNumber)]
        outputsDir = self.tempDir + os.sep + "outputs"
        outputsArg = ["-o", outputsDir, "-e", outputsDir]
        jobNameArg = ["-N", self.jobName]
        emailArg = []
        if self.emailAddress:
            emailArg = ["-M", self.emailAddress, "-m", self.emailConditions]
        holdArg = []
        if self.dependencies:
            dependencies = [str(dependency) for dependency in self.dependencies]
            holdArg = ["-hold_jid", ",".join(dependencies)]
        startingArgs = [genericRunnerPaths["qsub"], "-cwd"]
        schedulerCommandList = startingArgs + jobNameArg + limitsArgs + peArg + outputsArg + jobRangeArg + emailArg + holdArg
        commandToExecute = [genericRunnerPaths["python3"], genericRunnerPaths["arrayWrapper"], "-d", self.tempDir]
        commandToExecuteString = " ".join(commandToExecute)
        if self.mock:  #trap mock submissions here
            print("MOCK SUBMIT: " + " ".join(schedulerCommandList))
            print("MOCK COMMUNICATE: " + commandToExecuteString)
            print("Mock job, no number assigned")
            return 123456   #placeholder number
        successfullySubmitted = False
        retries = 0
        while not successfullySubmitted:
            child = subprocess.Popen(schedulerCommandList, stdout = subprocess.PIPE, stdin = subprocess.PIPE, stderr = subprocess.PIPE)
            childOut, childErr = child.communicate(input = commandToExecuteString.encode())
            childExitStatus = child.returncode
            qsubOut = childOut.decode().strip()
            qsubError = childErr.decode().strip()
            successfullySubmitted = childExitStatus == 0  #if it went through, it should exit status zero
            if successfullySubmitted:
                import re
                print("QSUB: %s" %(qsubOut))
                regex = re.search('Your job.* (\d+)', qsubOut)
                return int(regex.group(1))
            else:
                if retries >= maxRetries:
                    raise RuntimeError("Failed to qsub after %s attempts.\nQSUB OUT: %s\n QSUB ERR: %s" %(retries, qsubOut, qsubError))
                else:
                    print("Submission to qsub failed. Retrying.\nQSUB OUT: %s\n QSUB ERR: %s" %(qsubOut, qsubError))
                    retries += 1
        

        
    
        
    
    