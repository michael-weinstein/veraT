#!/usr/bin/env python3

import sys
import os
if not os.getcwd() + os.sep + "runners" in sys.path:
    sys.path.append(os.getcwd() + os.sep + "runners")
import programRunners
import genericRunners
import houseKeeping