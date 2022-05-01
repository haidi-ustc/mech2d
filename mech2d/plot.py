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

import os,string
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from monty.serialization import  loadfn,dumpfn
from mech2d.utils import prettyprint,box_center

debug=False
class Plot(object):

    def __init__(self,data='result.json',fmt='jpg',dpi=100):
        """
        Create an Plot obj. based on data

        Args: 
            data : Input data , should be a list or filename
        """
        self.fmt=fmt
        self.dpi=dpi
        if data=='result.json':
           self._data=self.load_result()      
        else:
           self._data=self.load_data(data)

    def load_result(self):
        """
        load json format data
        """
        with open('result.json','r') as f:
             d=json.loads(f.read())
        self.YV=d['YV']
        self.SS=d['SS']
        
    def load_data(self,data):
        """
        load json txt data or convert list to np.array
        """
        if os.path.isfile(data):
           return np.loadtxt(data)
        if isinstance(data,list):
           return np.array(data)

    def _set_fig(self,font_dict={},figsize=(8,6)):

        if font_dict:
           font = font_dict
        else:
           font =   {'family': 'serif', 'size': 14, 'weight':'bold'}
        plt.rc('font', **font)
        self.fig = plt.figure(figsize=figsize) # width , heighth
    
    def energy_strain_plot(self,font_dict={},figsize=(8,6),order=4,fname='Energy_Strain'):
        """
        plot energy strain figure 

        Args:
            font_dict (dict):   dict parameter for figure plot
            figsize (tuple):   dicide the figure size
            order (int):      fitting order for polynomial 
            fname (str):      output filename for figure
        Returns:
            None
        """
 
        self._set_fig(font_dict=font_dict,figsize=figsize)
        ax=self.fig.add_subplot(111)
        strain=self._data[:,0]*100
        energy=self._data[:,1]
        if debug:
           prettyprint(self._data)
        strain_min=strain.min()*1.01
        strain_max=strain.max()*1.01
        strain_vals=np.linspace(strain_min,strain_max,100)
        coeff=np.polyfit(strain.copy(),energy.copy(),order)
        poly=np.poly1d(coeff)
        energy_vals=poly(strain_vals)

        ax.scatter(strain,energy,marker='o',facecolor='r',edgecolor='k',alpha=0.7)
        ax.plot(strain_vals,energy_vals,linestyle='--',alpha=0.5)
        ax.set_xlabel('Lagrangian strain (%)') 
        ax.set_ylabel('Energy (eV)')
        plt.savefig(fname+'.'+self.fmt,format=self.fmt,dpi=self.dpi)
        #plt.show()
 
    def stress_strain_plot_nonefit(self,font_dict={},figsize=(16,16),fname='LStress_Strain'):
        """
        plot stress strain figure without fitting (for stress-strain curve calculation)

        Args:
            font_dict (dict):   dict parameter for figure plot
            figsize (tuple):   dicide the figure size
            order (int):      fitting order for polynomial 
            fname (str):      output filename for figure
        Returns:
            None
        """
        self._set_fig(font_dict=font_dict,figsize=figsize)
        if debug:
           prettyprint(self._data)
        strain=self._data[:,0]*100
        stress_XX=self._data[:,1]
        stress_YY=self._data[:,2]
        stress_XY=self._data[:,6]
        energy=self._data[:,7]
        strain_min=strain.min()*1.01
        strain_max=strain.max()*1.01
        Ys=[stress_XX,stress_YY,stress_XY,energy]
        labels=["XX",'YY',"XY"]
        loc=[[2,2,1],[2,2,2],[2,2,3],[2,2,4]]
        for ii,Y in enumerate(Ys):
            _loc=loc[ii]
            ax=self.fig.add_subplot(_loc[0],_loc[1],_loc[2])
            ax.plot(strain,Y,linestyle='-',marker='o',alpha=0.7)
            ax.set_xlabel('Lagrangian strain (%)') 
            if ii==3:
               ax.set_ylabel('Energy (eV)')
            else:
               ax.set_ylabel('%s Stress (N/m)'%(labels[ii]))
        plt.savefig(fname+'.'+self.fmt,format=self.fmt,dpi=self.dpi)

    def stress_strain_plot(self,font_dict={},figsize=(18,6),order=[3,3,3],fname='LStress_Strain'):
        """
        plot energy strain figure 

        Args:
            font_dict (dict):   dict parameter for figure plot
            figsize (tuple):   dicide the figure size
            order (int):      fitting order for polynomial 
            fname (str):      output filename for figure
        Returns:
            None
        """

        self._set_fig(font_dict=font_dict,figsize=figsize)
        if debug:
           prettyprint(self._data)
        strain=self._data[:,0]*100
        stress_XX=self._data[:,1]
        stress_YY=self._data[:,2]
        stress_XY=self._data[:,6]
        strain_min=strain.min()*1.01
        strain_max=strain.max()*1.01
        strain_vals=np.linspace(strain_min,strain_max,100)
        Ys=[stress_XX,stress_YY,stress_XY]
        labels=["XX",'YY',"XY"]
        for ii,Y in enumerate(Ys):
            ax=self.fig.add_subplot(1,3,ii+1)
            coeff=np.polyfit(strain.copy(),Y.copy(),order[ii])
            poly=np.poly1d(coeff)
            Y_vals=poly(strain_vals)
            ax.scatter(strain,Y,marker='o',facecolor='r',edgecolor='k',alpha=0.7)
            ax.plot(strain_vals,Y_vals,linestyle='--',alpha=0.5)
            ax.set_xlabel('Lagrangian strain (%)') 
            ax.set_ylabel('%s Stress (GPa)'%(labels[ii]))
        plt.savefig(fname+'.'+self.fmt,format=self.fmt,dpi=self.dpi)
        #plt.show()

    def polar_plot_EV(self,font_dict={},figsize=(8,6),skip=[1,1,1],fname='EV'):
        """
        plot the orientation dependent Young's Modulus and Poission's Ratio

        Args:
            font_dict (dict):   dict parameter for figure plot
            figsize (tuple):   dicide the figure size
            order (int):      fitting order for polynomial 
            fname (str):      output filename for figure
        Returns:
            None
        """
         
        self._set_fig(font_dict=font_dict,figsize=figsize)
        ax=self.fig.add_subplot(121, projection='polar')# ,facecolor="lightgoldenrodyellow")
          
        theta= self._data[0::skip[0],0]
        Y=self._data[0::skip[0],1]
        ax.plot(theta, Y, color="tab:orange", lw=1, ls="--", marker='h',alpha=0.6, label="$E$")
        ax.set_rlabel_position(10)  # get radial labels away from plotted line
        ax.tick_params(grid_color="palegoldenrod")
        angle = np.deg2rad(67.5)
        ax.legend(loc="lower left", bbox_to_anchor=(.5 + np.cos(angle)/2, .5 + np.sin(angle)/2))
        
        ax1=self.fig.add_subplot(122, projection='polar') #,facecolor="lightgoldenrodyellow")
        ax1.tick_params(grid_color="palegoldenrod")

        V=self._data[:,2]
        T=self._data[:,0]
        Vp_idx=np.where(V>0)
        Vn_idx=np.where(V<0)
        
        Vp=V[Vp_idx][0::skip[1]]
        theta_p=T[Vp_idx][0::skip[1]]
        
        Vn=V[Vn_idx][0::skip[2]]
        theta_n=T[Vn_idx][0::skip[2]]
 
        ax1.plot(theta_p, Vp, color="tab:green", lw=1, ls="--", marker='h',alpha=0.6, label=r"+ $\nu$")
        if len(Vn)>0:
           ax1.plot(theta_n, np.abs(Vn), color="tab:red", lw=0.1, marker='o',alpha=.5, label=r"- $\nu$")
        angle = np.deg2rad(67.5)
        ax1.legend(loc="lower left", bbox_to_anchor=(.5 + np.cos(angle)/2, .5 + np.sin(angle)/2))
        #ax1.set_rmax(2)
        ax1.set_rlabel_position(45)  # get radial labels away from plotted line
        #ax1.set_rticks([0.3, 0.60, 0.9])  # less radial ticks
        #ax.set_rticks([30,60,90])  # less radial ticks
        #ax.text(-0.25, 1.1, '('+string.ascii_lowercase[0]+')', transform=ax.transAxes,size=20, weight='bold')
        #ax1.text(-0.25, 1.1, '('+string.ascii_lowercase[1]+')', transform=ax1.transAxes,size=20, weight='bold')
        plt.tight_layout()
        plt.savefig(fname+'.'+self.fmt,format=self.fmt,dpi=self.dpi)
        
        #plt.show()

