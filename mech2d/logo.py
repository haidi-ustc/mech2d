#!/usr/bin/env python
#    mech2d is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.

#    mech2d is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public Licensea
#    along with mech2d.  If not, see <http://www.gnu.org/licenses/>.

from mech2d.utils import box_center
try:
   from mech2d._version import version
except:
   version="Unknow"
try:
   from mech2d._date import date
except:
   date="2021.09.05"

def logo():

    logo_list=     ['                     _     ___     _ ',
                    '                    | |   |__ \   | |',
                    ' _ __ ___   ___  ___| |__    ) |__| |',
                    "| '_ ` _ \ / _ \/ __| '_ \  / // _` |",
                    '| | | | | |  __/ (__| | | |/ /| (_| |',
                    '|_| |_| |_|\___|\___|_| |_|____\__,_|']

    box_center(ch='_',fill='_',sp=' ')
    box_center(ch=' ',fill=' ',sp='|')
    for lstr in logo_list:
        box_center(ch=lstr,fill=' ',sp='|')
    box_center(ch=' ',fill=' ',sp='|')
    box_center(ch=date+' '+version,fill=' ',sp='|')
    box_center(ch='_',fill='_',sp='|')
