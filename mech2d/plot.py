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


import numpy as np
import string
from monty.serialization import  loadfn,dumpfn


def Plot(object):
    def __init__(data='result.json'):
        if data=='result.json':
           self._data=self.load_result()      
        else:
           self._data=self.load_data(data)

    def load_result(self):
        with open('result.json','r') as f:
             d=json.loads(f.read())
        self.YV=d['YV']
        self.SS=d['SS']

    def load_data(data):
        if os.path.isfile(data):
           return np.loadtxt(data)
        if isinstance(data,list):
           return np.array(data)
        
    def polar_plot_EV(self,font_dict=None,figsize=None,skip=[1,1,1]):
         
        import matplotlib.pyplot as plt
        if font_dict:
           font = font_dict
        else:
           font =   {'family': 'serif', 'size': 14, 'weight':'bold'}
        plt.rc('font', **font)
       
        if figsize:
           pass
        else:
           figsize=(8,6) 
        fig = plt.figure(figsize=figsize) # width , heighth
        
        ax=fig.add_subplot(121, projection='polar')# ,facecolor="lightgoldenrodyellow")
        theta= f[0::skip[0],0]
        Y=f[0::skip[0],1]
        ax.plot(theta, Y, color="tab:orange", lw=1, ls="--", marker='h',alpha=0.6, label="$E$")
        ax.set_rlabel_position(10)  # get radial labels away from plotted line
        ax.tick_params(grid_color="palegoldenrod")
        angle = np.deg2rad(67.5)
        ax.legend(loc="lower left", bbox_to_anchor=(.5 + np.cos(angle)/2, .5 + np.sin(angle)/2))
        
        ax1=fig.add_subplot(122, projection='polar') #,facecolor="lightgoldenrodyellow")
        ax1.tick_params(grid_color="palegoldenrod")
        V=f[:,2]
        T=f[:,0]
        Vp_idx=np.where(V>0)
        Vn_idx=np.where(V<0)
        
        Vp=V[Vp_idx][0::skip[1]]
        theta_p=T[Vp_idx][0::skip[1]]
        
        Vn=V[Vn_idx][0::skip[2]]
        theta_n=T[Vn_idx][0::skip[2]]
 
        ax1.plot(theta_p, Vp, color="tab:green", lw=1, ls="--", marker='h',alpha=0.6, label=r"+ $\nu$")
        if len(Vn)>0:
           ax1.plot(theta_n, Vn, color="tab:red", lw=0.1, marker='o',alpha=.5, label=r"- $\nu$")
        angle = np.deg2rad(67.5)
        ax1.legend(loc="lower left", bbox_to_anchor=(.5 + np.cos(angle)/2, .5 + np.sin(angle)/2))
        #ax1.set_rmax(2)
        ax1.set_rlabel_position(45)  # get radial labels away from plotted line
        #ax1.set_rticks([0.3, 0.60, 0.9])  # less radial ticks
        #ax.set_rticks([30,60,90])  # less radial ticks
        ax.text(-0.25, 1.1, '('+string.ascii_lowercase[0]+')', transform=ax.transAxes,size=20, weight='bold')
        if len(Vn)>0:
           ax1.text(-0.25, 1.1, '('+string.ascii_lowercase[1]+')', transform=ax1.transAxes,size=20, weight='bold')
        plt.tight_layout()
        plt.savefig('EV.jpg',format='jpg',dpi=300)
        
        #plt.show()
def plot_elastic(args):
    
    plot=Plot()
