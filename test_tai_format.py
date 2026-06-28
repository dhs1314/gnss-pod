"""Test which UTC-TAI format Orekit v13 accepts."""
import os
os.environ['JAVA_HOME'] = r"C:\Program Files\JetBrains\PyCharm Community Edition 2024.3.5\jbr"

import orekit_jpype
orekit_jpype.initVM()

from java.io import File
from org.orekit.data import DataContext, DirectoryCrawler

data_dir = r"d:\prj\gnss_pod\data\orekit"

# Try different formats, test via DataProvidersManager
formats = [
    # Format 1: Standard IERS format with 6 columns
    ("IERS 6-col (YYYY MM DD Npl UT1-UTC TAI-UTC)",
     """ 1972  1  1  0  0 10
 1972  7  1  0  0 11
 1973  1  1  0  0 12
 1974  1  1  0  0 13
 1975  1  1  0  0 14
 1976  1  1  0  0 15
 1977  1  1  0  0 16
 1978  1  1  0  0 17
 1979  1  1  0  0 18
 1980  1  1  0  0 19
 1981  7  1  0  0 20
 1982  7  1  0  0 21
 1983  7  1  0  0 22
 1985  7  1  0  0 23
 1988  1  1  0  0 24
 1990  1  1  0  0 25
 1991  1  1  0  0 26
 1992  7  1  0  0 27
 1993  7  1  0  0 28
 1994  7  1  0  0 29
 1996  1  1  0  0 30
 1997  7  1  0  0 31
 1999  1  1  0  0 32
 2006  1  1  0  0 33
 2009  1  1  0  0 34
 2012  7  1  0  0 35
 2015  7  1  0  0 36
 2017  1  1  0  0 37
"""),
    # Format 2: Just YYYY MM DD offset (4 columns)
    ("4-col YYYY MM DD offset",
     """1972 1 1 10
1972 7 1 11
1973 1 1 12
1974 1 1 13
1975 1 1 14
1976 1 1 15
1977 1 1 16
1978 1 1 17
1979 1 1 18
1980 1 1 19
1981 7 1 20
1982 7 1 21
1983 7 1 22
1985 7 1 23
1988 1 1 24
1990 1 1 25
1991 1 1 26
1992 7 1 27
1993 7 1 28
1994 7 1 29
1996 1 1 30
1997 7 1 31
1999 1 1 32
2006 1 1 33
2009 1 1 34
2012 7 1 35
2015 7 1 36
2017 1 1 37
"""),
    # Format 3: MJD offset (2 columns)
    ("2-col MJD offset",
     """41317 10
41499 11
41683 12
42048 13
42413 14
42778 15
43144 16
43509 17
43874 18
44239 19
44786 20
45151 21
45516 22
46247 23
47161 24
47892 25
48257 26
48804 27
49169 28
49534 29
50083 30
50630 31
51179 32
53736 33
54832 34
56109 35
57204 36
57754 37
"""),
]

# Also test different filenames
filenames = ["UTC-TAI.history", "tai-utc.dat"]

for name, content in formats:
    for fname_base in filenames:
        fname = os.path.join(data_dir, fname_base)
        with open(fname, 'w') as f:
            f.write(content)

        DM = DataContext.getDefault().getDataProvidersManager()
        DM.clearProviders()
        DM.clearLoadedDataNames()
        DM.resetFiltersToDefault()

        data_file = File(data_dir)
        DM.addProvider(DirectoryCrawler(data_file))

        try:
            from org.orekit.frames import FramesFactory
            from org.orekit.utils import IERSConventions
            itrf = FramesFactory.getITRF(IERSConventions.IERS_2010, True)
            print(f"{name} + {fname_base}: OK")
        except Exception as e:
            msg = str(e)
            if "UTC-TAI" in msg:
                print(f"{name} + {fname_base}: UTC-TAI FAILED")
            elif "EOP" in msg:
                print(f"{name} + {fname_base}: EOP FAILED")
            else:
                print(f"{name} + {fname_base}: OTHER ERROR: {msg[:100]}")
