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


        
    
        
    
    