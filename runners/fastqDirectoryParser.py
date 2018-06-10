#!/usr/bin/env python3


class FastqDirectory(object):
    
    def __init__(self, directory, followlinks = True, depthLimit = -1):
        import os
        if not directory.endswith(os.sep):
            directory += os.sep
        if not os.path.isdir(directory):
            raise FastqDirectoryError("FASTQ directory %s was not found." %(directory))
        self.fileList = []
        for path, directories, files in os.walk(directory, followlinks = followlinks):
            if depthLimit > -1:
                currentDepth = (path + os.sep).replace(directory, "").count(os.sep)
                if currentDepth > depthLimit:
                    continue
            for file in files:
                if file.lower().endswith(".fastq.gz") or file.lower().endswith(".fastq"):
                    fileObject = FastqFile(os.sep.join([path, file]))
                    self.fileList.append(fileObject)
    
    def find(self, criteria, useLastResult = False):
        if type(criteria) == str:
            self.lastResult = self.quickL0Search(criteria, useLastResult)
        elif type(criteria) in (list, tuple):
            self.lastResult = self.listSearch(list(criteria), useLastResult)
        elif type(criteria) == dict:
            self.lastResult = self.dictSearch(criteria, useLastResult)
        elif type(criteria) == int:
            self.lastResult = self.intSearch(criteria, useLastResult)
        else:
            raise FastqDirectoryError("Invalid criteria type passed for search.\nType: %s\nValue: " %(type(criteria), criteria))
        return self.lastResult
        
    def quickL0Search(self, criteria, useLastResult = False):
        results = []
        if useLastResult:
            if hasattr(self, "lastResult"):
                searchList = self.lastResult
            else:
                raise FastqDirectoryError("Unable to refine nonexistent results (no last result set to refine).")
        else:
            searchList = self.fileList
        for file in self.fileList:
            if criteria in file.splitSampleName[0]:
                results.append(file)
        return results
    
    def dictSearch(self, criteria, useLastResult = False):
        rawCriteria = criteria
        criteria = {}
        for rawKey in list(rawCriteria.keys()):
            try:
                key = int(rawKey)
                criteria[key] = rawCriteria[rawKey]
            except (TypeError, ValueError):
                raise FastqDirectoryError("At least one non-integer key value in search dictionary.\nProblem value: %n\nFull dictionary: %s" %(rawKey, rawCriteria))
        results = []
        if useLastResult:
            if hasattr(self, "lastResult"):
                searchList = self.lastResult
            else:
                raise FastqDirectoryError("Unable to refine nonexistent results (no last result set to refine).")
        else:
            searchList = self.fileList
        for file in self.fileList:
            passed = True
            for key in list(criteria.keys()):
                try:
                    if not file.splitSampleName[key].lower() == criteria[key].lower():
                        passed = False
                        break
                except IndexError:
                    passed = False
                    break
            if passed:
                results.append(file)
        return results
    
    def listSearch(self, criteria, useLastResult = False):
        enum = enumerate(criteria)
        criteriaDict = {}
        for item in enum:
            criteriaDict[item[0]] = item[1]
        return self.dictSearch(criteriaDict, useLastResult)  #lazy or smart, you decide
    
    def intSearch(self, criteria, useLastResult = False):
        results = []
        if useLastResult:
            if hasattr(self, "lastResult"):
                searchList = self.lastResult
            else:
                raise FastqDirectoryError("Unable to refine nonexistent results (no last result set to refine).")
        else:
            searchList = self.fileList
        for file in self.fileList:
            if len(file.splitSampleName) >= criteria:
                results.append(file)
        return results
    
    def getAllFilesForSample(self, criteria, returnMatrix = False):
        #organizing by sampleName, then lane, then pairedEnd
        #print([str(item) for item in self.fileList])
        searchResults = self.find(criteria)
        self.collisionCheck(searchResults)
        sampleTree = {}
        for file in searchResults:
            if not file.sampleName in sampleTree:
                sampleTree[file.sampleName] = {}
            if not file.lane in sampleTree[file.sampleName]:
                sampleTree[file.sampleName][file.lane] = {}
            sampleTree[file.sampleName][file.lane][file.pairedEnd] = file
            if len(sampleTree[file.sampleName][file.lane]) > 2:
                raise FastqDirectoryError("Found more than two paired end files going into %s lane %s\n%s" %(file.sampleName, file.lane, sampleTree))
        fileMatrix = []
        sampleIndex = 0
        for sampleKey in list(sampleTree.keys()):
            fileMatrix.append([])
            for laneKey in list(sampleTree[sampleKey].keys()):
                fileMatrix[sampleIndex].append(sampleTree[sampleKey][laneKey])
            sampleIndex += 1
        if not returnMatrix:
            return sampleTree
        else:
            return fileMatrix
        
    def collisionCheck(self, fileList):
        fileInfoDict = {}
        for file in fileList:
            if not file.infoTuple in fileInfoDict:
                fileInfoDict[file.infoTuple] = file
            else:
                raise FastqDirectoryError("File naming collision detected between two files below:\n%s\n%s" %(fileInfoDict[file.infoTuple], file))
        return True        
    

class FastqFile(object):
    
    def __init__(self, file):
        if file == False:
            self.filePath = False
        else:
            import os
            if not os.path.isabs(file):
                file = os.path.abspath(file)
            self.filePath = file
            self.splitFileAndDirectory(file)
            if self.gzip:
                if not self.isValidGzipFile(file):
                    raise RuntimeError("File %s has the gzip file ending of '.gz' but does not appear to be gzip encoded. Please look into this and correct it." %file)
            self.infoTuple = (self.sampleName, self.lane, self.pairedEnd)
        
    def __str__(self):
        return self.filePath
        
    def splitFileAndDirectory(self, file):
        import os
        sepSplit = file.split(os.sep)
        self.directory = os.sep.join(sepSplit[:-1])
        self.fileName = sepSplit[-1]
        self.validateFileName(self.fileName)
        self.analyzeFileName(self.fileName)
        
    def validateFileName(self, fileName):
        fileNameDotList = fileName.lower().split(".")
        if fileNameDotList[-1] == "gz":
            self.gzip = True
            if not len(fileNameDotList) == 3:
                raise FastqDirectoryError("File named %s does not appear to conform to naming convention [name].fastq.gz" %(fileName))
            if not fileNameDotList [-2] == "fastq":
                raise FastqDirectoryError("Invalid file passed as FastqFile object: %s" %(fileName))
        else:
            self.gzip = False
            if not len(fileNameDotList) == 2:
                raise FastqDirectoryError("File named %s does not appear to conform to naming convention [name].fastq" %(fileName))
            if not fileNameDotList [-1] == "fastq":
                raise FastqDirectoryError("Invalid file passed as FastqFile object: %s" %(fileName))

    def isValidGzipFile(self, filename):
        import gzip
        file = gzip.open(filename)
        try:
            testline = file.readline()
            file.close()
        except OSError:
            file.close() 
            return False
        return True

        
    def analyzeFileName(self, fileName):
        fileName = fileName.split(".")[0]
        fileNameSplit = fileName.split("_")
        self.sampleName = "_".join(fileNameSplit[:-4])
        self.tags = "_".join(fileNameSplit[-4:])
        self.analyzeFileTags(self.tags)
        self.analyzeSampleName(self.sampleName)
    
    def analyzeFileTags(self, fileTags):
        rawSampleNumber, rawLaneNumber, rawPairedEnd, rawFileNumber = fileTags.split("_")
        try:
            self.sampleNumber = int(rawSampleNumber.replace("S",""))
        except ValueError:
            pass
        self.lane = int(rawLaneNumber.replace("L",""))
        self.pairedEnd = int(rawPairedEnd.replace("R",""))
        self.fileNumber = int(rawFileNumber)
        
    def analyzeSampleName(self, sampleName, forceDelimiter = False):
        usingFoundNonWordNonUnderscore = True
        if not forceDelimiter:
            import re
            nonWordNonUnderscore = re.findall("\W", sampleName)
            if len(nonWordNonUnderscore) == 0:
                if "_" in sampleName:
                    delimiter = "_"
                else:
                    self.splitSampleName = [sampleName] #nothing to split on here
                    return True
            elif len(nonWordNonUnderscore) > 1:
                print("WARNING: Sample name in %s has multiple non-word dividing characters. Unable to analyze sample name." %(self.fileName))
                self.splitSampleName = re.split("\W", sampleName)
                return True
            else: #only other possibility would be 1
                delimiter = nonWordNonUnderscore[0]
        else:
            delimiter = forceDelimiter
        self.splitSampleName = sampleName.split(delimiter)
        
class FastqDirectoryError(Exception):
    
    def __init__(self, message):
        self.message = message
        
    def __str__(self):
        return repr(self.message)
    
def cleanDirectory(directory):
    import re
    import os
    if not os.path.isdir(directory):
        raise FileNotFoundError("Unable to find a directory at %s" %directory)
    def breakName(name):
        return name.split(".")
    def fixName(name):
        return ".".join(name)
    for path, directories, files in os.walk(directory):
        print("%s\t%s\t%s" %(path, directories, files))
        for file in files:
            if not (file.endswith(".fastq") or file.endswith("fastq.gz")):
                continue
            oldName = file
            newName = breakName(oldName)
            newName = [re.sub("\W", "_", part) for part in newName]
            newName = fixName(newName)
            print("%s %s %s" %(oldName, newName, oldName == newName))
            if not newName == oldName:
                oldPath = path + os.sep + oldName
                newPath = path + os.sep + newName
                print("Changing %s to %s" %(oldPath, newPath))
                os.replace(oldPath, newPath)

if __name__ == '__main__':
    import sys
    import time
    args = sys.argv
    if not args[1] == "clean":
        raise RuntimeError("Only permitted action is \"clean\" while %s was passed" %args[1])
    directory = args[2]
    print("Starting cleanup in 10 seconds (press CTRL-C to abort)...", end = "", flush = True)
    time.sleep(10)
    print("STARTING")
    cleanDirectory(directory)
    print("DONE")