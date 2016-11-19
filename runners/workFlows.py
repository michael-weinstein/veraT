#!/usr/bin/env python3

import os

class PipelineStart(object):
    
    def __init__(self, sampleName, tempDir = False):
        if tempdir:
            import os
            if not os.path.isdir(tempDir):
                os.mkdir(tempDir)
            else:
                if os.path.isfile(tempDir + os.sep + "jobs.pkl"):
                    raise RuntimeError("A jobs.pkl already exists in the specified temporary directory.  Please clear the directory before overwriting with the new jobs.")
        else:
            import runnerSupport
            tempDir = runnerSupport.createTempDir(sampleName)
        self.tempDir = tempDir
        self.jobs = runnerSupport.JobList(self.tempDir)
        
    def end(self):
        self.jobs.close()
        
class AlignMemAddReadGroups(object):
    
    def __init__(self, refGenomeFasta, pe1, pe2 = False, bwaCores = 2, clobber = None, outputDir = ""):
        import programRunners
        import houseKeeping
        import picardRunners
        self.align = programRunners.BWAmem(sampleName, refGenomeFasta, pe1, pe2, cores = cores, clobber = None, outputDir = outputDir)
        self.bwaCommand = self.align.bwaCommand
        self.samFile = self.align.samOut
        self.makeBam = programRunners.ViewSAMtoBAM(sampleName, self.samFile, refGenomeFasta, clobber = self.align.clobber, outputDir = outputDir)
        self.viewSamCommand = self.makeBam.viewCommand
        self.bamFile = self.makeBam.bamOut
        self.deleteSam = houseKeeping.Delete(self.samFile)
        self.deleteSamCommand = self.deleteSam.deleteCommand
        self.addReadGroups = picardRunners.AddReadGroups(sampleName, self.bamFile, rgid = readGroupID, rglb = readGroupLibrary, rgpl = readGroupPlatform, rgpu = barcode, clobber = self.makeBam.clobber, outputDir = outputDir)
        self.addReadGroupsCommand = self.addReadGroups.picardCommand