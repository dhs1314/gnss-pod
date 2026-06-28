"""Quick test: orekit-jpype VM initialization."""
import os
os.environ['JAVA_HOME'] = r"C:\Program Files\JetBrains\PyCharm Community Edition 2024.3.5\jbr"

import orekit_jpype
print("orekit_jpype imported")

orekit_jpype.initVM()
print("VM initialized successfully")

# Test basic Orekit classes
from org.orekit.frames import FramesFactory
print("FramesFactory imported")

from org.orekit.time import TimeScalesFactory
print("TimeScalesFactory imported")

from org.orekit.utils import Constants
print("Orekit Earth mu:", Constants.EIGEN5C_EARTH_MU)
print("All OK!")
