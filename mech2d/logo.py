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

from mech2d.utils import box_center
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
    box_center(ch='version 1.0.0',fill=' ',sp='|')
    box_center(ch='_',fill='_',sp='|')
