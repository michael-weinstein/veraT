#!/usr/bin/env python3

__all__ = ["dataViewer",
           "getPolypeptides",
           "getRNASupportMPileup",
           "hlaReader",
           "makeOncotatorOutput",
           "mutectReader",
           "netMHCAnalyzer",
           "netMHCReader",
           "oncotatorReader",
           "peptideListMaker",
           "proteinSequenceFinder",
           "tandemVariantCombine",
           "variantCombine",
           "variantDataHandler",
           "variantSupport",
           "varScanPositionPuller",
           "varScanReader",
           "vcfReader"]

from . import dataViewer
from . import getPolypeptides
from . import getRNASupportMPileup
from . import hlaReader
from . import makeOncotatorOutput
from . import mutectReader
from . import netMHCAnalyzer
from . import netMHCReader
from . import oncotatorReader
from . import peptideListMaker
from . import proteinSequenceFinder
from . import tandemVariantCombine
from . import variantCombine
from . import variantDataHandler
from . import variantSupport
from . import varScanPositionPuller
from . import varScanReader
from . import vcfReader

#import sys
#import os
# if not os.getcwd() + os.sep + "runners" + os.sep + "variantReaders" in sys.path:
#     sys.path.append(os.getcwd() + os.sep + "runners" + os.sep + "variantReaders")
# import variantDataHandler
# import varScanReader
# import mutectReader
# import varScanPositionPuller
# import vcfReader