"""Find correct tai-utc.dat format for Orekit v13."""
import os
os.environ['JAVA_HOME'] = r"C:\Program Files\JetBrains\PyCharm Community Edition 2024.3.5\jbr"

import orekit_jpype
orekit_jpype.initVM()

# Look at USNO/IERS standard formats
# The standard IERS tai-utc.dat format has fixed columns:
# Columns: 1-4: year, 6-8: month, 10-12: day, 14-15: (blank), 17-27: TAI-UTC
# Example:
# " 1972  1  1      10.0"
# " 1972  7  1      11.0"

data_dir = r"d:\prj\gnss_pod\data\orekit"
fname = os.path.join(data_dir, "tai-utc.dat")

# Format 1: IERS fixed-width format (columns)
format1 = """ 1972  1  1      10.0
 1972  7  1      11.0
 1973  1  1      12.0
 1974  1  1      13.0
 1975  1  1      14.0
 1976  1  1      15.0
 1977  1  1      16.0
 1978  1  1      17.0
 1979  1  1      18.0
 1980  1  1      19.0
 1981  7  1      20.0
 1982  7  1      21.0
 1983  7  1      22.0
 1985  7  1      23.0
 1988  1  1      24.0
 1990  1  1      25.0
 1991  1  1      26.0
 1992  7  1      27.0
 1993  7  1      28.0
 1994  7  1      29.0
 1996  1  1      30.0
 1997  7  1      31.0
 1999  1  1      32.0
 2006  1  1      33.0
 2009  1  1      34.0
 2012  7  1      35.0
 2015  7  1      36.0
 2017  1  1      37.0
"""

# Format 3: Try including the pre-1972 fractional values in IERS format
format3 = """ 1961  1  1     1.422818
 1961  8  1     1.372818
 1962  1  1     1.845858
 1963 11  1     1.945858
 1964  1  1     3.240130
 1964  4  1     3.340130
 1964  9  1     3.440130
 1965  1  1     3.540130
 1965  3  1     3.640130
 1965  7  1     3.740130
 1965  9  1     3.840130
 1966  1  1     4.313170
 1968  2  1     4.213170
 1972  1  1     10.0
 1972  7  1     11.0
 1973  1  1     12.0
 1974  1  1     13.0
 1975  1  1     14.0
 1976  1  1     15.0
 1977  1  1     16.0
 1978  1  1     17.0
 1979  1  1     18.0
 1980  1  1     19.0
 1981  7  1     20.0
 1982  7  1     21.0
 1983  7  1     22.0
 1985  7  1     23.0
 1988  1  1     24.0
 1990  1  1     25.0
 1991  1  1     26.0
 1992  7  1     27.0
 1993  7  1     28.0
 1994  7  1     29.0
 1996  1  1     30.0
 1997  7  1     31.0
 1999  1  1     32.0
 2006  1  1     33.0
 2009  1  1     34.0
 2012  7  1     35.0
 2015  7  1     36.0
 2017  1  1     37.0
"""

formats = [
    ("IERS fixed-width (no pre-1972)", format1),
    ("IERS fixed-width (with pre-1972)", format3),
]

from java.io import File
from org.orekit.data import DataContext, DirectoryCrawler
from org.orekit.time import TimeScalesFactory
from org.orekit.frames import FramesFactory
from org.orekit.utils import IERSConventions

for name, content in formats:
    # Write file
    with open(fname, 'w') as f:
        f.write(content)
    print(f"\n--- Testing: {name} ---")
    print(f"File content:\n{content[:200]}...")

    # Reset
    DM = DataContext.getDefault().getDataProvidersManager()
    DM.clearProviders()
    DM.clearLoadedDataNames()
    DM.resetFiltersToDefault()
    DM.addProvider(DirectoryCrawler(File(data_dir)))

    try:
        utc = TimeScalesFactory.getUTC()
        last = utc.getLastKnownLeapSecond()
        print(f"UTC OK! Last leap second: {last}")
    except Exception as e:
        print(f"UTC FAILED: {e}")

    try:
        itrf = FramesFactory.getITRF(IERSConventions.IERS_2010, True)
        print(f"ITRF OK!")
    except Exception as e:
        print(f"ITRF FAILED: {e}")
