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

import os
import shutil
import numpy as np
from mech2d.constants import Len

def sepline(ch='-',sp='-'):
    r'''
    seperate the output by '-'
    '''
    print(ch.center(Len,sp))

def box_center(ch='',fill='-',sp="-"):
    r'''
    put the string at the center of |  |
    '''
    strs=ch.center(Len,fill)
    print(sp+strs[1:len(strs)-1:]+sp)


def create_path (path,back=False) :
    if  path[-1] != "/":
        path += '/'
    if os.path.isdir(path) :
        if back:
           dirname = os.path.dirname(path)
           counter = 0
           while True :
               bk_dirname = dirname + ".bk%03d" % counter
               if not os.path.isdir(bk_dirname) :
                   shutil.move (dirname, bk_dirname)
                   break
               counter += 1
           os.makedirs (path)
           return path
        else:
           return path
    os.makedirs (path)
    return path

def sortlist(lst1, lst2):
    temp = copy.copy(lst1)

    lst3 = []
    lst4 = []

    temp.sort()

    for i in range(len(lst1)):
        lst3.append(lst1[lst1.index(temp[i])])
        lst4.append(lst2[lst1.index(temp[i])])

    return lst3, lst4

def prettyprint(c,precision=3):
    #c=np.array(c)
    
    fmt="{:>10."+str(precision)+"f} "
    row = c.shape[0]
    col = c.shape[1]
    for i in range(row):
        for j in range(col):
            print(fmt.format(c[i, j]), end=" ")
            if j == (col - 1):
                print(" ")

if __name__=='__main__':
   create_path('test')
   create_path('test',back=True)
   
