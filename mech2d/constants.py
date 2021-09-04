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

r'''
 This module used to control the convert the unit
'''

# basic constant
_e     = 1.602176565e-19# elementary charge
Planck=6.6260755e-34    # Planck's constant, in Js

# Distance units
Bohr   = 5.291772086e-11# a.u. to meter
Ryd2eV = 13.605698066   # Ryd to eV
Bohr2Ang = 0.529177249  # Conversion of length from bohr to angstrom
Ang2Bohr = 1/Bohr2Ang
Ang2M= 1.0e-10 
M2Ang= 1.0/Ang2M 

# Energy units
Hartree2kCal = 627.5095 # Hartree to kcal/mol conversion
kCal2Hartree = 1/Hartree2kCal

eV2kCal = 23.061        # Conversion of energy in eV to energy in kcal/mol
kCal2eV = 1/eV2kCal

Hartree2Joule = 4.3597482e-18   # Hatree to Joule conversion factor
Joule2Hartree = 1/Hartree2Joule

eV2Hartree = Hartree2kCal/eV2kCal
Hartree2eV = 1/eV2Hartree

eV2Joule = 1.602e-19 
Joule2eV= 1/eV2Joule

Ry2eV =0.073498618 
eV2Py=1/Ry2eV

Ry2Hartree=0.50
Hartree2Ry=1/Ry2Hartree

#Pressure
Pa2GPa=1.0e-9
RdyToGPa  = (_e*Ryd2eV)/(1e9*Bohr**3)    # Ryd/[a.u.^3] to GPa
eVToGPa   = _e/(1e9*Ang2M**3)
eVToNpm   = _e/(Ang2M**2)   # eV/A**2 --> N/m
#others
eps_zero=0.0
Len=80
