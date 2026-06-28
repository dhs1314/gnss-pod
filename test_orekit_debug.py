"""Test Orekit without clearing default providers."""
import os
os.environ['JAVA_HOME'] = r"C:\Program Files\JetBrains\PyCharm Community Edition 2024.3.5\jbr"

import orekit_jpype
orekit_jpype.initVM()

from java.io import File
from org.orekit.data import DataContext, DirectoryCrawler

data_dir = r"d:\prj\gnss_pod\data\orekit"

# Option 1: Don't clear, just add
print("=== Option 1: Add to defaults ===")
DM = DataContext.getDefault().getDataProvidersManager()
DM.addProvider(DirectoryCrawler(File(data_dir)))
print(f"Added data dir: {data_dir}")

from org.orekit.time import TimeScalesFactory
try:
    utc = TimeScalesFactory.getUTC()
    print(f"UTC OK: last leap second = {utc.getLastKnownLeapSecond()}")
except Exception as e:
    print(f"UTC FAILED: {e}")

from org.orekit.frames import FramesFactory
from org.orekit.utils import IERSConventions
try:
    itrf = FramesFactory.getITRF(IERSConventions.IERS_2010, True)
    print(f"ITRF (simple=True) OK: {itrf.getName()}")
except Exception as e:
    print(f"ITRF FAILED: {e}")

# Option 2: Clear and re-add
print("\n=== Option 2: Clear + re-add ===")
DM.clearProviders()
DM.clearLoadedDataNames()
DM.resetFiltersToDefault()
DM.addProvider(DirectoryCrawler(File(data_dir)))
print("Cleared and re-added")

try:
    utc = TimeScalesFactory.getUTC()
    print(f"UTC OK: last leap second = {utc.getLastKnownLeapSecond()}")
except Exception as e:
    print(f"UTC FAILED: {e}")

try:
    itrf = FramesFactory.getITRF(IERSConventions.IERS_2010, True)
    print(f"ITRF (simple=True) OK: {itrf.getName()}")
except Exception as e:
    print(f"ITRF FAILED: {e}")
