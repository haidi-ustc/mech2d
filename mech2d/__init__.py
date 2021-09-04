#!/usr/bin/env python
#    elastic2d is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.

#    elastic2d is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public Licensea
#    along with elastic2d.  If not, see <http://www.gnu.org/licenses/>.

import logging
import os


ROOT_PATH=__path__[0]
NAME="mech2d"
SHORT_CMD="m2d"
log = logging.getLogger(__name__)
log.setLevel(logging.INFO)
logf = logging.FileHandler(os.getcwd()+os.sep+SHORT_CMD+'.log', delay=True)
logf_formatter=logging.Formatter('%(asctime)s - %(levelname)s : %(message)s')
logf.setFormatter(logf_formatter)
log.addHandler(logf)

__author__    = "Haidi Wang"
__copyright__ = "Copyright 2021"
__status__    = "Development"
try:
    from ._version import version as __version__
except ImportError:
    __version__ = 'unkown'
try:
    from ._date import date as __date__
except ImportError:
    __date__ = 'unkown'

def info():
    """
        Show basic information about """+NAME+""", its location and version.
    """
    logo()

    print('DeepModeling\n------------')
    print('Version: ' + __version__)
    print('Date:    ' + __date__)
    print('Path:    ' + ROOT_PATH)
    print('')
    print('Dependency')
    print('------------')
    for modui in ['numpy', 'dpdata', 'pymatgen', 'monty', 'ase', 'custodian' ]:
         try:
             mm = __import__(modui)
             print('%10s %10s   %s' % (modui, mm.__version__, mm.__path__[0]))
         except ImportError:
             print('%10s %10s Not Found' % (modui, ''))
    print()

    # reference
#    print("""Reference
#------------
#Please cite:
#mech2d: A tool for calculating elastic properties of two-dimensional materials
#from first principles, Computer Physics Communications, 2021, XXXX.
#------------
#""")

