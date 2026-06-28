"""Test Orekit with local IERS C04 data (no downloads needed)."""
import os, sys
os.environ['JAVA_HOME'] = r"C:\Program Files\JetBrains\PyCharm Community Edition 2024.3.5\jbr"

import orekit_jpype
orekit_jpype.initVM()

from java.io import File
from org.orekit.data import DataContext, DirectoryCrawler

# Use the existing IERS C04 file as Orekit data source
# Orekit's EopC04FilesLoader can read the IERS C04 14 format
data_dir = r"d:\prj\gnss_pod\data\orekit"
os.makedirs(data_dir, exist_ok=True)

# Copy the IERS C04 file to orekit data dir if not already there
iers_src = r"d:\prj\gnss_pod\data\IERS\eopc04_IAU2000.txt"
iers_dst = os.path.join(data_dir, "eopc04_IAU2000.txt")
if os.path.exists(iers_src) and not os.path.exists(iers_dst):
    import shutil
    shutil.copy2(iers_src, iers_dst)
    print(f"Copied {iers_src} -> {iers_dst}")

# Setup Orekit data provider pointing to the directory
DM = DataContext.getDefault().getDataProvidersManager()
DM.clearProviders()
DM.clearLoadedDataNames()
DM.resetFiltersToDefault()

data_file = File(data_dir)
if data_file.exists():
    DM.addProvider(DirectoryCrawler(data_file))
    print(f"Data provider set to: {data_dir}")
else:
    print(f"ERROR: Data directory {data_dir} does not exist")
    sys.exit(1)

# Try to load EOP via Orekit
try:
    from org.orekit.frames import FramesFactory, LazyLoadedFrames
    from org.orekit.utils import IERSConventions

    # Try getting ITRF with IERS 2010
    itrf = FramesFactory.getITRF(IERSConventions.IERS_2010, True)
    print(f"ITRF (IERS 2010): {itrf}")
    print(f"ITRF name: {itrf.getName()}")

    # Try getting GCRF
    gcrf = FramesFactory.getGCRF()
    print(f"GCRF: {gcrf.getName()}")

    # Try time scales
    from org.orekit.time import TimeScalesFactory
    utc = TimeScalesFactory.getUTC()
    tt = TimeScalesFactory.getTT()
    print(f"UTC: {utc}")
    print(f"TT: {tt}")

    # Check leap seconds
    last_leap = utc.getLastKnownLeapSecond()
    print(f"Last known leap second: {last_leap}")

    # Build a test date
    from org.orekit.time import AbsoluteDate
    date = AbsoluteDate(2024, 4, 29, 12, 0, 0.0, utc)
    print(f"Test date: {date}")

    # Try frame transform
    from org.orekit.utils import PVCoordinates
    from org.hipparchus.geometry.euclidean.threed import Vector3D
    pos = Vector3D(7000000.0, 0.0, 0.0)  # ~7000 km
    vel = Vector3D(0.0, 7500.0, 0.0)     # ~7.5 km/s
    pv = PVCoordinates(pos, vel)
    print(f"Test PV in GCRF: pos={pos}, vel={vel}")

    # Transform GCRF -> ITRF
    transform = gcrf.getTransformTo(itrf, date)
    pv_itrf = transform.transformPVCoordinates(pv)
    print(f"Transformed PV in ITRF: pos={pv_itrf.getPosition()}, vel={pv_itrf.getVelocity()}")

    print("\n=== Orekit basic functionality OK! ===")

except Exception as e:
    import traceback
    print(f"ERROR: {e}")
    traceback.print_exc()
