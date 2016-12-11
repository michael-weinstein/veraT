#!/usr/bin/env python3

global genericRunnerPaths
houseKeepingRunnerPaths = {}

class Delete(object):
    
    def __init__(self, file):
        if type(file) == list or type(file) == tuple:
            file = " ".join(file)
        self.deleteCommand = "rm -f " + file
        
class Move(object):
    
    def __init__(self, source, destination):
        self.moveCommand = "mv " + source + " " + destination
        
class Gzip(object):
    def __init__(self, file):
        self.gzipCommand = "gzip " + file
        
class Capstone(object):
    def __init__(self, sampleName, outputDir = ""):
        if not outputDir:
            import os
            outputDir = os.getcwd() + os.sep
        capstoneCommand = "echo \"SAMPLE " + sampleName + " DONE\" > " + outputDir + sampleName + ".capstone ; sleep 180"


        
    
        
    
    