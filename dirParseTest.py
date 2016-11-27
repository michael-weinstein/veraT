#!/usr/bin/env python3

import runners.fastqDirectoryParser

dirParse = runners.fastqDirectoryParser.FastqDirectory("sampleData")
sampleMatrix = dirParse.getAllFilesForSample("JMS")
sampleTree = dirParse.getAllFilesForSample("JMS", returnTree = True)
hashSearch = dirParse.find({2:"022213"})
listSearch = dirParse.find(["JMS","normal"])
intSearch = dirParse.find(3)
print("something")