#!/usr/bin/env python3

import os
global genericRunnerPaths
houseKeepingRunnerPaths = {}
global thisProgramPath
thisProgramPath = (os.path.abspath(__file__))

class CheckArgs(object):
    
    def __init__(self):
        import argparse
        import os
        parser = argparse.ArgumentParser()
        parser.add_argument("-m", "--mode", help = "Mode to run in")
        parser.add_argument("-f", "--files", help = "Files to act upon")
        rawArgs = parser.parse_args()
        mode = rawArgs.mode
        if not mode:
            raise RuntimeError("No mode specified, no job to do.")
        else:
            self.mode = mode.upper()
        files = rawArgs.files
        if "," in files:
            files = files.split(",")
        else:
            files = [files]
        self.files = files

class Delete(object):  #slightly fancier than a simple rm command, will try to remove associated bai files when removing a bam file
    
    def __init__(self, file):
        if type(file) == list or type(file) == tuple:
            file = ",".join(file)
        self.deleteCommand = "/u/local/apps/python/3.4.3/bin/python3 %s -m delete -f %s" %(thisProgramPath, file)
        
class Move(object):
    
    def __init__(self, source, destination):
        self.moveCommand = "mv " + source + " " + destination
        
class Gzip(object):
    def __init__(self, file):
        self.gzipCommand = "gzip " + file
        
class Capstone(object):
    def __init__(self, sampleName, outputDir = ""):
        import os
        if not outputDir:
            outputDir = os.getcwd() + os.sep
        self.capstoneCommand = "echo \"SAMPLE " + sampleName + " DONE\" > " + outputDir + os.sep + sampleName + ".capstone ; sleep 180"

def fileDeleter(files):
    if type(files) == str:
        if "," in files:
            files = files.split(",")
        else:
            files = [files]
    import os
    for file in files:
        if os.path.isfile(file):
            os.remove(file)
        if file.endswith(".bam"):
            if os.path.isfile(file + ".bai"):
                os.remove(file + ".bai")
            elif os.path.isfile(file[:-3] + "bai"):
                os.remove(file[:-3] + "bai")
    
if __name__ == '__main__':
    args = CheckArgs()
    if args.mode == "DELETE":
        fileDeleter(args.files)
    else:
        raise RuntimeError("Invalid mode set: %s" %(args.mode))
    quit()
        
    
    