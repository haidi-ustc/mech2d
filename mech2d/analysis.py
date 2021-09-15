#!/usr/bin/env python

#    mech2d is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public V3 License as published by
#    the Free Software Foundation.

#    mech2d is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with mech2d.  If not, see <http://www.gnu.org/licenses/>.

import numpy as np
from numpy import sin,cos,pi
from mech2d.constants import Len
from mech2d.utils import prettyprint,sepline,box_center
from mech2d.bravais import Bravais2D
from mech2d.plot import Plot
from pymatgen.core import Structure

class Analysis(object):
    """
     A module to analysis the Young's modulus , Poisson's ratio and Shear modulus
    """
    def __init__(self, structure, elastic_tensor,plot=False,approach=None):

        self.C2d = elastic_tensor
        self.structure = structure 
        self.plot = plot
        self.approach=approach if approach else ''

    def get_brav_lattice(self):
        return Bravais2D(self.structure)

    @property
    def lattice_type(self):
        brav = self.get_brav_lattice()
        if brav.lattice_type == "Oblique":
           return 'O'
        elif brav.lattice_type == "CenteredRectangular":
           return 'CR'
        elif brav.lattice_type == "Rectangular":
           return 'R'
        elif brav.lattice_type == "Square":
           return 'S'
        elif brav.lattice_type == "Hexagonal":
           return 'H'
        else:
           raise RuntimeError('Make sure you have a standard conventional\n cell with vacuum layer along z direction.')

    def get_gamma(self):
        """
        return the layer modulus.
        """
        return 0.25 * (self.C2d[0][0] + self.C2d[1][1] + 2 * self.C2d[0][1])

    def get_Y10(self):
        """
        return 2D Young's modulus or in-plane stiffness: Y[10] = [c11c22 - c12^2]/[c22]
        """
        return (self.C2d[0][0] * self.C2d[1][1] - self.C2d[0][1] * self.C2d[0][1]) / self.C2d[1][1]

    def get_Y01(self):
        """
        2D Young's modulus or in-plane stiffness: Y[01] = [c11c22 - c12^2]/[c11]
        """
        return (self.C2d[0][0] * self.C2d[1][1] - self.C2d[0][1] * self.C2d[0][1]) / self.C2d[0][0]

    def get_nu10(self):
        """
        2D Poisson's ratio;  nu10 = c12/c22
        """
        return self.C2d[0][1] / self.C2d[1][1]

    def get_nu01(self):
        """
        2D Poisson's ratio; nu01 = c12/c11
        """
        return self.C2d[0][1] / self.C2d[0][0]

    def get_G2d(self):
        """
        2D shear modulus; G2d = C66

        Args:
            None
        Returns:
            shear modulus 
        """
        return self.C2d[5][5]

    def get_stability(self):
        if self.lattice_type == 'O':
            if self.C2d[0,0]>0 and self.C2d[0,0]*(self.C2d[1,1])>self.C2d[0,1]**2 and np.linalg.det(self.C2d)>0:
               return 'Stable'
            else:
               return 'Unstable'
        elif self.lattice_type == 'CR' or self.lattice_type == 'R':
            if self.C2d[0,0]>0 and self.C2d[0,0]*(self.C2d[1,1])>self.C2d[0,1]**2 and self.C2d[5,5]>0:
               return 'Stable'
            else:
               return 'Unstable'
        elif self.lattice_type == 'S':
            if self.C2d[0,0]>0 and self.C2d[0,0]>abs(self.C2d[0,1]) and self.C2d[5,5]>0:
               return 'Stable'
            else:
               return 'Unstable'
        elif self.lattice_type == 'H':
            if self.C2d[0,0]>0 and self.C2d[0,0]>abs(self.C2d[0,1]):
               return 'Stable'
            else:
               return 'Unstable'
        else:
           raise RuntimeError('ERROR: Unknown 2D bravais lattice')

    def references(self):
        box_center(ch="References",fill='-',sp='-')
        refs = "Phys. Rev. B 85, 125428 (2012)\n"
        refs+= "Acta Mechanica 223, 2591-2596 (2012)\n"
        refs+= "Comput. Mater. Sci. 68, 320 (2013)\n"
        refs+= "Mech. Mater. 64, 135 (2013)\n"
        refs+= "2D Mater. 6, 048001 (2019)"
        print(refs)

    def summary(self):
        box_center(ch="Elastic properties summary",fill='-',sp='-')
        print('')
        box_center(ch="Stiffness Tensor of %s system (N/m)"%(self.structure.formula.replace(" ","")),fill='-',sp='-')
        prettyprint(self.C2d)
        tmp=np.zeros((6,6))
        tmp[2,2]=100
        tmp[3,3]=100
        tmp[4,4]=100
        self.S2d = np.linalg.inv(self.C2d+tmp)
        if self.get_stability()=='Stable':
           box_center(ch="Compliance Tensor (m/N)",fill='-',sp='-')
           self.S2d[2,2]=0
           self.S2d[3,3]=0
           self.S2d[4,4]=0
           prettyprint(self.S2d,precision=5)
           sepline()
           print('Writing orientation-dependent E and v ...')
           self.get_EV()
        sepline()
        print("2D layer modulus (N/m)         :   %10.3f " % self.get_gamma())
        print("2D Young's modulus Y[10] (N/m) :   %10.3f " % self.get_Y10())
        print("2D Young's modulus Y[01] (N/m) :   %10.3f " % self.get_Y01())
        print("2D Shear modulus G (N/m)       :   %10.3f " % self.get_G2d())
        print("2D Poisson ratio v[10]         :   %10.3f " % self.get_nu10())
        print("2D Poisson ratio v[01]         :   %10.3f " % self.get_nu01())
        print("Stability                      :   %10s   " % self.get_stability())
        self.references()
        box_center(ch='-',fill='-',sp='-')

    @property
    def v_zz(self):
        """
        return vzz=C12/C22
        """
        return self.C2d[0][1]/self.C2d[1][1]

    @property
    def d1(self):
        """
        return d1=C11/C22+1-(C11*C22-C12**2)/C22/C66;
        """
        return self.C2d[0][0]/self.C2d[1][1]+1 - \
               (self.C2d[0][0]*self.C2d[1][1]-self.C2d[0][1]**2)/ \
               self.C2d[1][1]/self.C2d[5][5]

    @property
    def d2(self):
        """
        return d2=-(2*C12/C22-(C11*C22-C12**2)/C22/C66);
        """
        return -1*(2*self.C2d[0][1]/self.C2d[1][1]-\
                (self.C2d[0][0]*self.C2d[1][1]-self.C2d[0][1]**2)/ \
                 self.C2d[1][1]/self.C2d[5][5])

    @property
    def d3(self):
        """
        return d3 =C11/C22
        """
        return self.C2d[0][0]/self.C2d[1][1]

    @property
    def Y_zz(self):
        """
        return Y_zz = C11*C22-C12**2)/C22
        """
        return (self.C2d[0][0]*self.C2d[1][1]-self.C2d[0][1]**2)/self.C2d[1][1]

    def get_EV(self,npoints=360,fname='EV_theta.dat'):
        theta=np.linspace(0,2*pi,360)
        E=self.Y_zz/((cos(theta))**4+self.d2*(cos(theta))**2.*(sin(theta))**2+self.d3*(sin(theta))**4)
        V=(self.v_zz*(cos(theta))**4-self.d1*(cos(theta))**2.*(sin(theta))**2+self.v_zz*(sin(theta))**4)/ \
          ((cos(theta))**4+self.d2*(cos(theta))**2.*(sin(theta))**2+self.d3*(sin(theta))**4);
        res=np.vstack((theta,E,V)).T
        np.savetxt(fname,res,fmt='%10.6f %10.6f %10.6f')
        if self.plot:
           try:
              _plot=Plot(data=fname)
              _plot.polar_plot_EV(fname=self.approach+'-EV.jpg')
           except:
              print('WARNING: Plot failed, skip !!!')

if __name__=='__main__':
   st=Structure.from_file('POSCAR')
   c=np.random.randn(6,6)*10
   c=c*c.T
   A=Analysis(st,c)
   A.references()
   A.summary()
