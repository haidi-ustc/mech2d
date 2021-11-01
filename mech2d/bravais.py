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

import math
import numpy as np
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


class Bravais2D:
    """ 
    2D bravais lattice obj to check the lattice type
 
    Args
        a: (float) The magnitude of the first primitive vector (default is 1.0).
        b: (float) The magnitude of the second primitive vector (default is 1.0).
        angle: (float) The angle between the two primitive vectors; can't be 0 or 180 degrees (default is 120.0).
        centered: (bool) True if the lattice is a centered rectangular (default is False).
        numpoints: (int) The number of desired points to plot and must be a square number larger than 4;
        will be the number of 'non-centered' points if centered rectangular lattice (default is 16).
        plot: (bool) If True, will plot the lattice (default is True).
        a_vec: (numpy array) The first primitive vector.
        b_vec: (numpy array) The second primitive vector.
        lattice: (str) The name of the type of Bravais lattice (depends on a, b, angle, and centered).
        unit_cell_area: (float) The area of the unit cell.
    Functions:
        plot: Creates a 2D scatter plot.
    """

    def __init__(self, struct , eps_r=0.2, eps_a=1.0, numpoints=16):
        """
        :param a: (float) The magnitude of the first primitive vector (default is 1.0).
        :param b: (float) The magnitude of the second primitive vector (default is 1.0).
        :param angle: (float) The angle between the two primitive vectors; can't be 0 or 180 degrees (default is 120.0).
        :param degrees: (bool) If True, the angles are in degrees and if False, the angles are in radians
        (default is True).
        :param centered: (bool) True if the lattice is a centered rectangular (default is False).
        :param numpoints: (int) The number of desired points to plot and must be a square number larger than 4;
        will be the number of 'non-centered' points if centered rectangular lattice (default is 16).
        :param plot: (bool) If True, will plot the lattice (default is True).
        """
        self.struct= struct
        lattice= struct.lattice.abc
        self.a = lattice[0]
        self.b = lattice[1]
        self.eps_r= eps_r
        self.eps_a= eps_a
        self.centered = self._get_centered()
        self._numpoints = numpoints
        self.angle = struct.lattice.gamma

    def _get_centered(self):
        #spa=SpacegroupAnalyzer(self.struct)
        if len(self.struct) != len(self.struct.get_primitive_structure()):
            centered=True
        else:
            centered=False
        return centered

    @property
    def a_vec(self):
        return np.array([self.a, 0])

    @property
    def b_vec(self):
        return np.array([self.b*math.cos(self.angle), self.b*math.sin(self.angle)])

    @property
    def lattice_type(self):
        #print(self.centered)
        #print(self.struct.lattice.parameters)
        #if (abs(self.angle - 90)>  self.eps_a) and (not self.centered):
        if (abs(self.a - self.b)> self.eps_r) and (abs(self.angle - 90)>  self.eps_a) and (not self.centered):
            return "Oblique"
        elif (abs(self.a - self.b)> self.eps_r) and (abs(self.angle - 90)<=  self.eps_a) and self.centered:
            return "CenteredRectangular"
        elif (abs(self.a - self.b)> self.eps_r) and (abs(self.angle - 90)<=  self.eps_a) and (not self.centered):
            return "Rectangular"
        elif (abs(self.a - self.b)<= self.eps_r) and ((abs(self.angle - 120)<= self.eps_a) or (abs(self.angle - 60) <= self.eps_a ) ) and (not self.centered):
            return "Hexagonal"
        #elif (abs(self.a - self.b)<= self.eps_r) and (abs(self.angle - 90)<=self.eps_a) and (not self.centered):
        elif (abs(self.a - self.b)<= self.eps_r) and (abs(self.angle - 90)<=self.eps_a) :
            return "Square"
        elif (abs(self.angle - 90)>  self.eps_a) and (not self.centered):
            return "Oblique"
        else:
            raise Exception("Invalid combination of a, b, angle, and/or centered.")

    @property
    def unit_cell_area(self):
        return self.a*self.b*math.sin(self.angle)

    @property
    def numpoints(self):
        val = round(self._numpoints**0.5, 0)**2
        if val == self._numpoints and val >= 9:
            return self._numpoints
        else:
            raise Exception("numpoints must be a square number greater than or equal to 9.")

    def __find_points(self):
        """ Finds all the x and y coordinates of the lattice.
        :return: (list(list)) x, y
        """

        def f(start, stop, x_list, y_list):
            for j in np.arange(start, stop):
                for i in np.arange(start, stop):
                    vec = i*self.a_vec + j*self.b_vec
                    x_list.append(vec[0])
                    y_list.append(vec[1])
            return x_list, y_list

        p = int(self.numpoints**0.5)
        x, y = f(0, p, [], [])
        if self.centered:
            x, y = f(0.5, p - 1, x, y)
        return x, y

    def __unit_cell(self, x, y):
        """ Finds the x and y coordinates for the unit cell.
        :param x: (list) The x coordinates of the lattice points.
        :param y: (list) The y coordinates of the lattice points.
        :return: (list(list()) x, y
        """

        root = int(self.numpoints**0.5)
        return (x[0], x[1], x[1+root], x[root], x[0]), (y[0], y[1], y[1+root], y[root], y[0])

    def plot(self):
        import matplotlib.pyplot as plt
        """ Creates a 2D scatter plot. """

        x, y = self.__find_points()
        fig = plt.figure()
        ax = fig.add_subplot(111)
        angle = "%.2f"%(self.angle)
        title = "Bravais Lattice: " + self.lattice_type
        variables = "\n|a| = " + "%.3f"%(self.a) + ", |b| = " + "%.3f"%(self.b) + ", \u03b8  =  " + angle + "\u00b0"
        scaling = "\n(Axes may be scaled differently)"
        ax.set_title(title + variables + scaling)
        ax.scatter(x, y, label="Lattice Points")
        ax.plot(*self.__unit_cell(x, y), color="darkorange", label="Unit Cell")
        plt.legend(loc="best")#.set_draggable(True)
        plt.show()

if __name__=='__main__':
   from pymatgen.core import Structure
   from glob import glob
   fs=glob("../tests/configs/*.vasp")
   structs=[Structure.from_file(f) for f in fs]
   for i,struct in enumerate(structs[:-2]):
       print('-'* 20)
       print(fs[i])
       brav2 = Bravais2D(struct)
       #print(struct.lattice)
       print(brav2.lattice_type)
   brav2.plot()
