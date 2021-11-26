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

import json
import os,sys
import numpy as np
from monty.json import MSONable
from pymatgen.core import Structure
from monty.serialization import loadfn,dumpfn

from mech2d.bravais import Bravais2D
from mech2d.utils import create_path
from mech2d.calculation.vasp import VASP
from mech2d import log
from mech2d.constants import Len, eVToNpm
from mech2d.utils import prettyprint,box_center
from mech2d.analysis import Analysis
from mech2d.plot import Plot
from mech2d.logo import logo

#--------------------------------------------------------------------------------------------------
# Define the global variables
# Lagrangian strain dict.
Ls_Dic={                       
'01':[ 1., 0., 0., 0., 0., 0.],
'02':[ 1., 1., 0., 0., 0., 0.],
'03':[ 0., 1., 0., 0., 0., 0.],
'04':[ 0., 0., 0., 0., 0., 2.],
'05':[ 1., 0., 0., 0., 0., 2.],
'06':[ 0., 1., 0., 0., 0., 2.]}

Ls_str={                                     
'01':'(  eta,  0.0,  0.0,  0.0,  0.0,  0.0)',
'02':'(  eta,  eta,  0.0,  0.0,  0.0,  0.0)',
'03':'(  0.0,  eta,  0.0,  0.0,  0.0,  0.0)',
'04':'(  0.0,  0.0,  0.0,  0.0,  0.0, 2eta)',
'05':'(  eta,  0.0,  0.0,  0.0,  0.0, 2eta)',
'06':'(  0.0,  eta,  0.0,  0.0,  0.0, 2eta)'}

Direct_Dict={
'xx': '01',
'bi': '02',
'yy': '03',
'xy': '04'
}

LT_Dic = {              
'O'  :'Oblique',
'R'  :'Rectangular',
'CR' :'rectangular',
'H'  :'Hexagonal',
'S'  :'Square'
} 

def_fmt1="%d"
def_fmt2="%03d"
felastic='Mech2D.json'
fresult='result.json'

#--------------------------------------------------------------------------------------------------
class Elastic(MSONable):

    """
    A elastic object for initializaion and post-process.  
    """

    def __init__(self,structure,approach='energy',properties='elc', workdir=None,verbose=False):
        """
        Create an Elastic obj. from structure with given condition
        Args: 
            structure : Input structure file name or Pymatgen Structure obj.
            approach (str): Two approaches are avaliable for elastic constant calculation, 'stress' or 'energy'  
            properties (str): What kind of properties do you want to calculate ? elastic constant or  stress-strain curve, supported
                        values are : 'elc' or 'ssc'
            workdir (str): The working directory, default is None
            verbose (bool): Output the verbose information for debug
        """
        if isinstance(structure,str) and os.path.isfile(structure):
           self._struct = Structure.from_file(structure)
        elif isinstance(structure,dict):
           self._struct = Structure.from_dict(structure)
        elif isinstance(structure,Structure):
           self._struct = structure
        else:
           raise RuntimeError('structure must be file name or Pymatgen Structure obj')
        
        self._approach = approach
        self._properties = properties
        self._verbose = verbose
        if workdir is None:
           self.workdir =os.path.join(os.getcwd(), self.properties+'_'+self.approach)
        else:
           self.workdir=workdir

    def __str__(self):
        ret = '-'*Len+'\n'
        if self.verbose:
           ret+= repr(self.structure)
        else:
           ret+= repr(self.structure.lattice)
    
        ret+= '\n'+'-'*Len
        ret+= "\napproach    : %s "%self.approach
        ret+= "\nproperties  : %s "%self.properties
        ret+= "\nlattice     : %s "%LT_Dic[self.lattice_type]
        ret+= "\nnumber SOEC : %s "%self.number_elastic_tensor
        ret+= '\n'+'-'*Len
        return ret
 
    def __repr__(self):
        return str(self)

    @property
    def verbose(self):
        """
        Verbose information
        """
        return self._verbose

    @verbose.setter
    def verbose(self,val):
        self._verbose=val

    @property
    def properties(self):
        """
        Property for calculation
        """
        return self._properties
      
    @property
    def approach(self):
        """
        Approach for calculation
        """
        return self._approach

    @property
    def structure(self):
        """
        Structure object 
        """
        return self._struct

    def get_brav_lattice(self):
        return Bravais2D(self.structure)
    
    @property
    def lattice_type(self):
        """
        Two dimensional lattice type
        """
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
           raise RuntimeError('Make sure you have a standard conventional\n cell with vacuum layer along z direction')

    @property
    def C_mult_eta_Matrix(self):
        """
        The relation between lagrangian stress and strain : tao=C*yeta
        """
        if self.lattice_type == 'O':
           #'c11 c12 c22 C16 C26 c66'
           return np.mat(
                        [
                        [1,0,0,0,0,0],
                        [0,1,0,0,0,0],
                        [0,0,0,1,0,0],
                        [1,1,0,0,0,0],
                        [0,1,1,0,0,0],
                        [0,0,0,1,1,0],
                        [1,0,0,2,0,0],
                        [0,1,0,0,2,0],
                        [0,0,0,1,0,2],
                        ]
                        )
        elif self.lattice_type == 'CR' or self.lattice_type == 'R':
           #'c11 c12 c22 c66'
           return np.mat([
                          [1,1,0,0],
                          [0,1,1,0],
                          [0,0,0,0],
                          [1,0,0,0],
                          [0,1,0,0],
                          [0,0,0,2]
                         ]
                          )
        elif self.lattice_type == 'S':
           # 'c11 c12 c66'
           return np.mat([[1,0,0],
                          [0,1,0],
                          [0,0,2]]
                          )
        elif self.lattice_type == 'H':
           # 'c11 c12'
           return np.mat([[1,0],
                          [0,1],
                          [1,-1]]
                          )
        else: 
           raise RuntimeError('ERROR: Unknown 2D bravais lattice')
    @property
    def number_elastic_tensor(self):
        """
        number of the independent elastic tensor 
        """
        if self.lattice_type == 'O':
           return 6
        elif self.lattice_type == 'CR' or self.lattice_type == 'R':
           return 4
        elif self.lattice_type == 'S':
           return 3
        elif self.lattice_type == 'H':
           return 2
        else: 
           raise RuntimeError('ERROR: Unknown 2D bravais lattice')

    @property
    def lagrangian_strain_list(self):
       """
       lagrangian strain list
       """
       if  self.approach == 'energy':

           if self.lattice_type == 'O':
               Lag_strain_list = ['01','02','03','04','05','06']
           elif self.lattice_type == 'CR' or self.lattice_type == 'R':
               Lag_strain_list = ['01','02','03','04']
           elif self.lattice_type == 'S':
               Lag_strain_list = ['01','02','04']
           elif self.lattice_type == 'H':
               Lag_strain_list = ['01','02']
           else:
              raise RuntimeError('ERROR: Unknown 2D bravais lattice')
    
       elif  self.approach =='stress':

           if self.lattice_type == 'O':
               Lag_strain_list = ['01','02','05']
           elif self.lattice_type == 'CR' or self.lattice_type == 'R':
               Lag_strain_list = ['02','05']
           elif self.lattice_type == 'S':
               Lag_strain_list = ['05']
           elif self.lattice_type == 'H':
               Lag_strain_list = ['05']
           else: 
              raise RuntimeError('ERROR: Unknown 2D bravais lattice')
       else:
           raise RuntimeError('ERROR: Unknown approach for elastic constant calcuator')
       return Lag_strain_list
          
    def as_dict(self) :
        """
        Json-serialization dict representation of the Elastic.

        """
        d = {
            "@module": self.__class__.__module__,
            "@class": self.__class__.__name__,
            "structure": self.structure,
            "approach": self.approach,
            "properties": self.properties,
            "lattice": self.lattice_type,
            "SOEC": self.number_elastic_tensor,
            "verbose":self.verbose,
            "workdir":self.workdir
        }
        return d
   
    @classmethod
    def from_dict(cls, d):
        """
        Reconstitute a Molecule object from a dict representation created using
        as_dict().

        Args:
            d (dict): dict representation of Elastic

        Returns:
            Elastic object
        """
        return cls(d['structure'],approach=d['approach'], properties=d['properties'],
                   workdir=d['workdir'], verbose = d['verbose'])

    def to(self,filename,option={},fmt=None,indent=4):

        """
        Outputs the Elastic to a file.

        Args:
            fmt (str): Format to output to. Defaults to JSON
            filename (str): Output filename
            option (dict): The optional parameters in dict will be append to Elastic.as_dict()
            indent (int): Used for format output
        Returns:
            None
        """

        d=self.as_dict()
        d.update(option) 
        if fmt and filename.endswith(fmt):
           filename+='.'+fmt
        dumpfn(d,filename,indent=indent) 
   
    def calc_stress_strain(self,skip=False,finput='input.yaml',plot=False,fmt='jpg',dpi=100):
        """
        Calculate the stress train curve.

        Args:
            skip (str): skip the data parsing by Mech2d parser
            finput (str): input file name for basic parameter about machine, task, resources, and code
            plot (bool): whether plot the figure, default is False
        """


        if self.approach=="energy":
           raise RuntimeError("ERROR: Stress Strain curve calculation only support stress approach")

        pwd=os.getcwd()
        inputs=loadfn(os.path.join(pwd,finput))
        code=inputs['code']['name'].upper()
        with open(os.path.join(self.workdir,felastic),'r') as f:
             elat2d=json.loads(f.read())
        numb_points=elat2d['numb_points']

        if 'max_lag_strain' in elat2d.keys():
            max_lag_strain=elat2d['max_lag_strain']
            lag_strain_range=None
        else:
            max_lag_strain=None
            lag_strain_range=elat2d['lag_strain_range']
            
        workdir=elat2d["workdir"]
        V0=elat2d["V0"]
        direction=elat2d['direction']
        height=self.structure.lattice.c 
        _cell=self.structure.lattice.matrix.copy()
        area = np.linalg.norm(np.cross(_cell[0],_cell[1]))

        for i in direction:
            def_key  = 'Def_'+i
            Ls_list= Ls_Dic[Direct_Dict[i]]
            Defn = os.path.join(pwd, self.workdir, def_key)
            sigs=[]
            taos=[]
            for j in range(1,numb_points+1):
                Defn_num = os.path.join(Defn, def_key+'_'+def_fmt2%j)
                if (os.path.exists(Defn_num)):

                    if code=='VASP': 
                       CONV=-1*height/100 # for VASP
                       ret0 = VASP.get_stress(Defn_num)
                       ret1 = VASP.get_energy(Defn_num)
                       if ret0 and ret1:
                           sig = np.array(ret0)
                           energy = float(ret1)
                       else:
                           print("ERROR: cannot parse energy or stress from:\n%s"%Defn_num)
                           os._exit(0)

                    sig[2][2]=0
                    sig[0][2]=0
                    sig[2][0]=0
                    sig[2][1]=0
                    sig[1][2]=0

                    if max_lag_strain:
                       rs=np.linspace(-max_lag_strain,max_lag_strain,numb_points)
                    else:
                       rs=np.linspace(lag_strain_range[0],lag_strain_range[1],numb_points)
                    r=rs[j-1]

                    le = np.zeros(6)
                    for k in range(6):
                        le[k] = Ls_list[k]

                    Lv = r*le
                    eps_matrix=self._get_eps_matrix(Lv.copy())
                    def_matrix = np.eye(3) + eps_matrix

                    # calculate the Lagrangian stress from phyical stress and physical strain
                    dm  = def_matrix.copy()
                    idm = np.linalg.inv(dm)
                    tao = np.linalg.det(dm)*np.dot(idm,np.dot(sig,idm))
                    tao=tao*CONV
                    sig=sig*CONV
                    sigs.append([r,sig[0][0],sig[1][1],sig[2][2],sig[1][2],sig[0][2],sig[0][1],energy])
                    taos.append([r,tao[0][0],tao[1][1],tao[2][2],tao[1][2],tao[0][2],tao[0][1],energy])
                else:
                    raise RuntimeError("ERROR: Directory %s not exits"%Defn_num)
            fL=os.path.join(Defn,def_key+'_Lagrangian_Stress.dat')
            fP=os.path.join(Defn,def_key+'_Physical_Stress.dat')
            if self.verbose:
               box_center("Physical Stress")
               prettyprint(np.array(sigs))
               box_center("Lagrangian Stress")
               prettyprint(np.array(taos))
            _fmt='%14.8f  '*7+'%16.8f'
            _header="#Phy. strain          XX           YY           ZZ           YZ           XZ           XY (N/m)  energy (eV)" 
            np.savetxt(fP, np.array(taos), delimiter=" ", header=_header, fmt=_fmt, comments='')
            _header="#Lag. strain          XX           YY           ZZ           YZ           XZ           XY (N/m)  energy (eV)" 
            np.savetxt(fL, np.array(sigs), delimiter=" ", header=_header, fmt=_fmt, comments='')
            if plot: 
               try:
                  _plot=Plot(data=fL,fmt=fmt,dpi=dpi)
                  _plot.stress_strain_plot_nonefit(fname=os.path.join(Defn,os.path.basename(Defn)+'_Lagrangian_Stress'))
               except:
                  print('WARNING: Plot failed, skip !!!')
                
#--------------------------------------------------------------------------------------------------
    def calc_elastic_constant(self,poly_order=4,skip=False,finput='input.yaml',plot=False,fmt='jpg',dpi=100):

        """
        Calculate the elstic constant using the different approach.

        Args:
            poly_order (int): the order of polynomial for fitting
            skip (str): skip the data parsing by Mech2d parser
            finput (str): input file name for basic parameter about machine, task, resources, and code
            plot (bool): whether plot the figure, default is False
        """

        if self.approach == 'energy':
           if isinstance(poly_order,int):
              poly_order=[poly_order]*len(self.lagrangian_strain_list)
           elif isinstance(poly_order,list):
              if len(poly_order) !=len(self.lagrangian_strain_list):
                 poly_order=[poly_order[0]]*len(self.lagrangian_strain_list)
           else:
              raise RuntimeError('ERROR: Polynomial order must be integer or list of integer with length = %d '
                                 %len(len(self.lagrangian_strain_list)))
        elif self.approach=='stress':
           if isinstance(poly_order,int):
              poly_order=[poly_order]*6
           elif isinstance(poly_order,list):
              if len(poly_order) !=6:
                 poly_order=[poly_order[0]]*6
           else:
              raise RuntimeError('ERROR: Polynomial order must be integer or list of integer with length = %d '
                                %len(len(self.lagrangian_strain_list)))
             
        else:
           raise RuntimeError('ERROR: Unkonwn approach! The supported are "energy" or "stress"')

        pwd=os.getcwd()
        inputs=loadfn(os.path.join(pwd,finput))
        code=inputs['code']['name'].upper()
        with open(os.path.join(self.workdir,felastic),'r') as f:
             elat2d=json.loads(f.read())
        numb_points=elat2d['numb_points']
        max_lag_strain=elat2d['max_lag_strain']
        workdir=elat2d["workdir"]
        V0=elat2d["V0"]
        
        #print(poly_order)
        height=self.structure.lattice.c 
        _cell=self.structure.lattice.matrix.copy()
        area = np.linalg.norm(np.cross(_cell[0],_cell[1]))
        
        #print(height,V0/area)
        if self.approach=='energy':
           C=self.calc_elastic_constant_from_energy(poly_order=poly_order,skip=skip,
                                                    numb_points=numb_points,
                                                    max_lag_strain=max_lag_strain,
                                                    workdir=workdir,
                                                    code=code,
                                                    plot=plot,fmt=fmt,dpi=dpi)
           C = eVToNpm*C/area # from eV/a^2 to N/m
        elif self.approach=='stress':
           C=self.calc_elastic_constant_from_stress(poly_order=poly_order,skip=skip,
                                                    numb_points=numb_points,
                                                    max_lag_strain=max_lag_strain,
                                                    workdir=workdir,
                                                    code=code,
                                                    plot=plot,fmt=fmt,dpi=dpi)

           if self.verbose:
              print('height %.3f'%height)
           C = -C*height/100 # from eV/a^2 to N/m
        else:
           raise RuntimeError('ERROR: Unkonwn approach! The supported are "energy" or "stress"')

        # unit convert 
        ana = Analysis(self.structure,C,plot=plot,approach=self.approach)
        ana.summary(fmt=fmt,dpi=dpi)

#--------------------------------------------------------------------------------------------------
    def calc_elastic_constant_from_stress(self,poly_order,skip,numb_points,max_lag_strain,workdir,code,plot,fmt='jpg',dpi=100):
        """
        Calculate the elstic constant based on stress-strain approach.
        """
        for i in range(1,len(self.lagrangian_strain_list)+1):
            Ls_list= Ls_Dic[self.lagrangian_strain_list[i-1]]
            Defn  = os.path.join(workdir,'Def_'+def_fmt1%i)
            if  not os.path.exists(Defn):
                raise RuntimeError('ERROR: The '+ Defn +' directory not found.')
            #os.chdir(Defn)
            fL=os.path.join(Defn,os.path.basename(Defn)+'_Lagrangian_Stress.dat')
            fP=os.path.join(Defn,os.path.basename(Defn)+'_Physical_Stress.dat')
            if skip: 
               # Sometimes user can parsing the output file by themselves 
               print("Skip parsing: \n%s \n%s"%(fL,fP) )
            else:
               with open(fL, 'w') as fidL, open(fP, 'w') as fidP:
                    fidL.write('#Lag. strain          XX           YY           ZZ           YZ           XZ           XY \n')
                    fidP.write('#Phy. strain          XX           YY           ZZ           YZ           XZ           XY \n')
                    for j in range(1, numb_points+1):
                        Defn_num = os.path.join(Defn, os.path.basename(Defn)+'_'+def_fmt2%j)

                        if (os.path.exists(Defn_num)):

                            if code=='VASP': 
                               ret = VASP.get_stress(Defn_num)
                               if ret:
                                   sig = ret
                               else:
                                   print("ERROR: cannot parse energy from:\n%s"%Defn_num)
                                   os._exit(0)
                              
                            sig[2][2]=0
                            sig[0][2]=0
                            sig[2][0]=0
                            sig[2][1]=0
                            sig[1][2]=0
                            s = j-(numb_points+1)/2
                            r = 2*max_lag_strain*s/(numb_points-1)

                            le = np.zeros(6)
                            for k in range(6):
                                le[k] = Ls_list[k]

                            Lv = r*le
                            eps_matrix=self._get_eps_matrix(Lv.copy())
                            def_matrix = np.eye(3) + eps_matrix

                            # calculate the Lagrangian stress from phyical stress and physical strain
                            dm  = def_matrix.copy()
                            idm = np.linalg.inv(dm)
                            tao = np.linalg.det(dm)*np.dot(idm,np.dot(sig,idm))

                            if (s==0): r=0.0001

                            if (r>0):
                                strain ='+%12.10f'%r
                            else:
                                strain = '%13.10f'%r

                            fidP.write(strain +'   '+'%10.8f'%sig[0][0]\
                                              +'   '+'%14.8f'%sig[1][1]\
                                              +'   '+'%14.8f'%sig[2][2]\
                                              +'   '+'%14.8f'%sig[1][2]\
                                              +'   '+'%14.8f'%sig[0][2]\
                                              +'   '+'%14.8f'%sig[0][1]+'\n')

                            fidL.write(strain +'   '+'%14.8f'%tao[0][0]\
                                              +'   '+'%14.8f'%tao[1][1]\
                                              +'   '+'%14.8f'%tao[2][2]\
                                              +'   '+'%14.8f'%tao[1][2]\
                                              +'   '+'%14.8f'%tao[0][2]\
                                              +'   '+'%14.8f'%tao[0][1]+'\n')
        
                        else:
                            raise RuntimeError("ERROR: Directory %s not exits"%Defn_num)
        if self.verbose:
           box_center('poly_order')
           print(poly_order)

        sigma=[] 
        for i in range(1,len(self.lagrangian_strain_list)+1):
            Defn  = os.path.join(workdir,'Def_'+def_fmt1%i)
            fL=os.path.join(Defn,os.path.basename(Defn)+'_Lagrangian_Stress.dat')
            if plot: 
               try:
                  _plot=Plot(data=fL,fmt=fmt,dpi=dpi)
                  _plot.stress_strain_plot(order=poly_order,
                                           fname=os.path.join(Defn,os.path.basename(Defn)+'_Lagrangian_Stress'))
               except:
                  print('WARNING: Plot failed, skip !!!')
            ret=np.loadtxt(fL,skiprows=1) 
            ld="strain          XX           YY           ZZ           YZ           XZ           XY".split()
            for idx in range(6):
                if idx in [2,3,4]:
                   continue
                if self.verbose:
                   box_center("fitting data: %s & %s "%(ld[0],ld[idx+1]))
                   prettyprint(np.vstack((ret[:,0],ret[:,idx+1])),precision=5)
                coeff = np.polyfit(ret[:,0], ret[:,idx+1], poly_order[idx])
                sigma.append(float(coeff[poly_order[idx]-1]))

        if self.verbose:
           box_center(ch='sigma')
           print(sigma)

        Matrix=self.C_mult_eta_Matrix
       
        sigma = np.array(sigma)
        
        if self.verbose:
           box_center(' Matrix')
           prettyprint(Matrix,precision=4)
        _C     = np.linalg.lstsq(Matrix,sigma,rcond=None)[0]

        if self.verbose:
           box_center(ch='_C')
           print(_C)
    
        C     = np.zeros((6,6))
        if self.lattice_type == 'O':
             C[0,0]=_C[0]
             C[0,1]=_C[1]
             C[1,0]=_C[1]
             C[1,1]=_C[2]
             c[0,5]=_C[3]
             c[5,0]=_C[3]
             c[1,5]=_C[4]
             c[5,1]=_C[4]
             C[5,5]=_C[5]

        elif self.lattice_type == 'CR' or self.lattice_type == 'R':
             C[0,0]=_C[0]
             C[0,1]=_C[1]
             C[1,0]=_C[1]
             C[1,1]=_C[2]
             C[5,5]=_C[3]
           
        elif self.lattice_type == 'S':
             C[0,0]=_C[0]
             C[0,1]=_C[1]
             C[1,0]=_C[1]
             C[1,1]=_C[0]
             C[5,5]=_C[2]

        elif self.lattice_type == 'H':
             C[0,0]=_C[0]
             C[0,1]=_C[1]
             C[1,0]=_C[1]
             C[1,1]=_C[0]
             C[5,5]=0.5*(_C[0]-_C[1])
            
        else: 
           raise RuntimeError('ERROR: Unknown 2D bravais lattice.')
        return C

#--------------------------------------------------------------------------------------------------
    def calc_elastic_constant_from_energy(self,poly_order,skip,numb_points,max_lag_strain,workdir,code,plot,fmt='jpg',dpi=100):
        """
        Calculate the elstic constant based on energy-strain approach.
        """

        for i in range(1,self.number_elastic_tensor+1):
            Defn  = os.path.join(workdir,'Def_'+def_fmt1%i)
            if  not os.path.exists(Defn):
                raise RuntimeError('ERROR: The '+ Defn +' directory not found.')
            # use the text format, easy use for users with other DFT code
            fEV=os.path.join(Defn,os.path.basename(Defn)+'_Energy.dat')

            if skip: 
               # Sometimes user can parsing the output file by themselves 
               print("Skip parsing: \n%s"%fEV)
            else:
               with open(fEV, 'w') as fid:
                    for j in range(1, numb_points+1):

                        Defn_num = os.path.join(Defn, os.path.basename(Defn)+'_'+def_fmt2%j)

                        if (os.path.exists(Defn_num)):
                            if code=='VASP': 
                               ret=VASP.get_energy(Defn_num)
                               if ret:
                                   energy = "%12.10f"%(float(ret))
                               else:
                                   print("ERROR: cannot parse energy from:\n%s"%Defn_num)
                                   os._exit(0)

                            s = j-(numb_points+1)/2
                            r = 2*max_lag_strain*s/(numb_points-1)

                            if (s==0): r=0.0001

                            if (r>0):
                                strain ='+%12.10f'%r
                            else:
                                strain = '%13.10f'%r
                            fid.write(strain+'   '+energy+'\n')

        A2 = []
        for i in range(1,self.number_elastic_tensor+1):
            Defn  = os.path.join(workdir,'Def_'+def_fmt1%i)
            fEV=os.path.join(Defn,os.path.basename(Defn)+'_Energy.dat')
            if plot: 
               try:
                  _plot=Plot(data=fEV,fmt=fmt,dpi=dpi)
                  _plot.energy_strain_plot(order=poly_order[i-1],
                                           fname=os.path.join(Defn,os.path.basename(Defn)+'_Energy_Strain'))
               except:
                  print('WARNING: Plot failed, skip !!!')
            ret=np.loadtxt(fEV) 
            if self.verbose:
               box_center("fitting data @ %s "%os.path.basename(Defn))
               prettyprint(np.vstack((ret[:,0],ret[:,1]-ret[int(ret.shape[0]/2),1])).T,precision=5)
            coeffs = np.polyfit(ret[:,0], ret[:,1],poly_order[i-1])
            # shift or not are both fine
            #coeffs = np.polyfit(ret[:,0], ret[:,1]-ret[int(ret.shape[0]/2),1], poly_order[i-1])  
            A2.append(coeffs[poly_order[i-1]-2])            

        if len(A2) != self.number_elastic_tensor:
           raise RuntimeError('ERROR: The number of *_Energy.dat is NOT equal to ' + str(self.number_elastic_tensor))

        C = np.zeros((6,6))
        if self.lattice_type == 'O':
            C[0,0]=2*A2[0]
            C[1,1]=2*A2[2]
            C[0,1]=A2[1]-A2[0]-A2[2]
            C[1,0]=C[0,1]
            C[5,5]=0.5*A2[3]
            C[0,5]=0.5*(A2[4]-A2[0]-A2[3])
            C[5,0]=C[0,5]
            C[1,5]=0.5*(A2[4]-A2[0]-A2[2])
            C[5,1]=C[1,5]

        elif self.lattice_type == 'CR' or self.lattice_type == 'R':
            C[0,0]=2*A2[0]
            C[1,1]=2*A2[2]
            C[0,1]=A2[1]-A2[0]-A2[2]
            C[1,0]=C[0,1]
            C[5,5]=0.5*A2[3]

        elif self.lattice_type == 'S':
            C[0,0]=2*(A2[0])
            C[0,1]=A2[1]-2*A2[0]
            C[1,0]=C[0,1]
            C[1,1]=C[0,0]
            C[5,5]=0.5*A2[2]
            
        elif self.lattice_type == 'H':
            C[0,0]=2*(A2[0])
            C[0,1]=A2[1]-2*A2[0]
            C[1,0]=C[0,1]
            C[1,1]=C[0,0]
            C[5,5]=0.5*(C[0,0]-C[0,1])

        else: 
           raise RuntimeError('ERROR: Unknown 2D bravais lattice')
        return C

#--------------------------------------------------------------------------------------------------
    @staticmethod
    def _get_eps_matrix(Lv):
        """
        Calculate the physical strain according to the Lagrangian strain based on iteration method.

        Args:
           Lv (array): Lagrangian strain array
        Returns:
           (array): physical strain array
        """
        # eps= eta - 0.5 * eps  , this matrix equation can be solved by iteration method
        eta_matrix      = np.zeros((3,3))
        eta_matrix[0,0] = Lv[0]
        eta_matrix[0,1] = Lv[5]/2.
        eta_matrix[0,2] = Lv[4]/2.

        eta_matrix[1,0] = Lv[5]/2.
        eta_matrix[1,1] = Lv[1]
        eta_matrix[1,2] = Lv[3]/2.

        eta_matrix[2,0] = Lv[4]/2.
        eta_matrix[2,1] = Lv[3]/2.
        eta_matrix[2,2] = Lv[2]

        eps_matrix = eta_matrix

        norm       = 1.0
        while( norm > 1.e-10 ):
            x          = eta_matrix - np.dot(eps_matrix, eps_matrix)/2.
            norm       = np.linalg.norm(x - eps_matrix)
            eps_matrix = x
        return eps_matrix

#--------------------------------------------------------------------------------------------------
    def set_deformed_structure_for_ss(self,numb_points,max_lag_strain=None,lag_strain_range=None,direction=['xx'],back=True):

        """
        set the deformed structure for stress-strain curve calculation

        Args:
            numb_points (int): Number of deformed structures
            max_lag_strain (float): The maximum Lagrangian strain, the range will be set to [-max_lag_strain,+max_lag_strain]
            lag_strain_range (array): Set the Lagrangian strain manually
            direction (list): Which direction to calculate the stress-strain curve
            back (bool) : back the old folder. Default True 
        Returns:
            None
        """

        if (numb_points < 2):
           raise RuntimeError('ERROR: Too few deformed structure, increase numb_points')

        if self.approach == 'energy':
           raise RuntimeError("ERROR: Stress Strain curve calculation only support stress approach")

        if self.approach == 'stress':
           interval = 0.00001

        if max_lag_strain:
           if (numb_points%2 == 0):
               numb_points   += 1
               print( 'The number of the Deformed structures change to  '+ str(numb_points) )
           else:
               print( 'The number of the Deformed structures is  '+ str(numb_points))
           ptn = int((numb_points-1)/2)
           if (max_lag_strain/ptn <= interval):
              raise RuntimeError('ERROR: The interval of the strain values is < '+ str(interval) +\
                  ' Choose a larger maximum Lagrangian strain'\
                  'or a less number of deformed structures.\n')
           loop_list = range(-ptn, ptn+1)
           rs=np.linspace(-max_lag_strain,max_lag_strain,numb_points)
        else:
           ptn = numb_points
           rs=np.linspace(lag_strain_range[0],lag_strain_range[1],ptn)
           if (abs(rs[1]-rs[0]) <= interval):
              raise RuntimeError('ERROR: The interval of the strain range is < '+ str(interval) +\
                  ' Choose a larger Lagrangian strain range'\
                  'or a less number of deformed structures.\n')
           loop_list = range(1,ptn+1)

        if self.verbose:
           box_center("ratio list")
           print(len(rs),len(loop_list))
           print(loop_list)   
        cell  = self._struct.lattice.matrix.copy()
        V0   = abs(np.linalg.det(cell))
        pwd = os.getcwd() 

        cell_old= cell.copy()
        create_path(self.workdir,back=back)
        def_params={}

        cont1= 0
        for i in direction:
            Ls_list= Ls_Dic[Direct_Dict[i]]
            _Ls_str= Ls_str[Direct_Dict[i]]
            if self.verbose:
               box_center('Lagrangian strain List')
               print(Ls_list)
               print(_Ls_str)
            cont1  = cont1 + 1
            def_key  = 'Def_'+i
            Defn = os.path.join(pwd, self.workdir,'Def_'+i)

            create_path(Defn,back=back)
            cont2 = 0
               
            def_structs_list = []
            for s,r in zip(loop_list,rs):

                def_structs_dict = {'Path':None,'eta':None,'Cell':None}
                
                Ls = np.zeros(6)
                for k in range(6):
                    Ls[k] = Ls_list[k]
                Lv = r*Ls
                 
                if self.verbose:
                   print(Lv)

                Lv_matrix=np.zeros((3,3))
                Lv_matrix[0][0]=Lv[0]
                Lv_matrix[1][1]=Lv[1]
                Lv_matrix[2][2]=Lv[2]
                Lv_matrix[2][1]=Lv[3]/2
                Lv_matrix[1][2]=Lv[3]/2
                Lv_matrix[2][0]=Lv[4]/2
                Lv_matrix[0][2]=Lv[4]/2
                Lv_matrix[0][1]=Lv[5]/2
                Lv_matrix[1][0]=Lv[5]/2

                eps_matrix = self._get_eps_matrix(Lv.copy())

                #--- Calculating the cell_new matrix ---------------------------------------------------------
                i_matrix   = np.eye(3)
                def_matrix = i_matrix + eps_matrix
                cell_new      = np.dot(cell_old, def_matrix)

                if self.verbose:
                   print('raitio %.3f    Lv_matrix        eta_mtrix          def_matrix'%r)
                   #box_center(ch='raitio %.3f'%r,fill=' ',sp=' ')
                   prettyprint(np.hstack((Lv_matrix,eps_matrix,def_matrix)),precision=5)
                #------------------------------------------------------------------------------------------
                cont2 = cont2 + 1
                Defn_cont2  = os.path.join(Defn, Defn.split('/')[-1]+'_'+def_fmt2%cont2)
        
                def_structs_dict['Path'] = Defn_cont2
                def_structs_dict['eta']  = r
                def_structs_dict['Cell'] = cell_new.tolist()
                def_structs_list.append(def_structs_dict)

                #--- Writing the structure file -----------------------------------------------------------
                create_path(Defn_cont2,back=False)
                new_struct=Structure(cell_new,self.structure.species,self.structure.frac_coords)
                new_struct.to("POSCAR",os.path.join(Defn_cont2,Defn_cont2.split('/')[-1] +'.vasp'))
            def_params[def_key]={'Path':Defn,'LagrangianStrain':Ls_str[Direct_Dict[i]],'DeformedCell':def_structs_list}
        #--------------------------------------------------------------------------------------------------
        dumpfn(def_params,os.path.join(self.workdir,'Deformed_Parameters.json'),indent=4)

        if max_lag_strain:
           self.to(os.path.join(self.workdir,felastic),
                   option={"numb_points":numb_points,"max_lag_strain":max_lag_strain,
                           "direction":direction,"V0":V0})
        else:
           self.to(os.path.join(self.workdir,felastic),
                   option={"numb_points":numb_points,"lag_strain_range":lag_strain_range,
                           "direction":direction,"V0":V0})

    def set_deformed_structure_for_elc(self,numb_points,max_lag_strain,back=True):
        """
        set the deformed structure for elastic constant calculation by using 'stress' or 'energy' approach

        Args:
            numb_points (int): Number of deformed structures
            max_lag_strain (float): The maximum Lagrangian strain, the range will be set to [-max_lag_strain,+max_lag_strain]
            back (bool) : back the old folder. Default True 
        """

        if (numb_points < 5):
           raise RuntimeError('ERROR: Too few deformed structure, increase numb_points')
        if (numb_points > 100):
           raise RuntimeError('ERROR: Too much deformed structure, decrease numb_points')
    
        if (numb_points%2 == 0):
            numb_points   += 1
            print( 'The number of the Deformed structures change to  '+ str(numb_points) )
        else:
            print( 'The number of the Deformed structures is  '+ str(numb_points))
    
        ptn = int((numb_points-1)/2)

        if self.approach == 'energy':
           interval = 0.0001
        if self.approach == 'stress':
           interval = 0.00001
        if (max_lag_strain/ptn <= interval):
           raise RuntimeError('ERROR: The interval of the strain values is < '+ str(interval) +\
               ' Choose a larger maximum Lagrangian strain'\
               'or a less number of deformed structures.\n')

        cell  = self._struct.lattice.matrix.copy()
        V0   = abs(np.linalg.det(cell))

        pwd = os.getcwd() 
    
        cell_old= cell.copy()
        create_path(self.workdir,back=back)
        def_params={}

        cont1= 0
        for i in self.lagrangian_strain_list:
            Ls_list= Ls_Dic[i]
            _Ls_str= Ls_str[i]
            if self.verbose:
               box_center('Lagrangian strain List')
               print(Ls_list)
               print(_Ls_str)
            cont1  = cont1 + 1
            def_key  = 'Def_'+def_fmt1%cont1
            Defn = os.path.join(pwd, self.workdir,'Def_'+def_fmt1%cont1)

            create_path(Defn,back=back)
        
            cont2 = 0
            
            loop_list = range(-ptn, ptn+1)
               
            def_structs_list = []
            for s in loop_list:

                def_structs_dict = {'Path':None,'eta':None,'Cell':None}
                r = max_lag_strain*s/ptn
                
                if (s==0):
                    if (self.approach == 'Energy'): r = 0.0001
                    if (self.approach == 'Stress'): r = 0.00001
        
                Ls = np.zeros(6)
                for k in range(6):
                    Ls[k] = Ls_list[k]
                Lv = r*Ls
                 
                if self.verbose:
                   print(Lv)

                Lv_matrix=np.zeros((3,3))
                Lv_matrix[0][0]=Lv[0]
                Lv_matrix[1][1]=Lv[1]
                Lv_matrix[2][2]=Lv[2]
                Lv_matrix[2][1]=Lv[3]/2
                Lv_matrix[1][2]=Lv[3]/2
                Lv_matrix[2][0]=Lv[4]/2
                Lv_matrix[0][2]=Lv[4]/2
                Lv_matrix[0][1]=Lv[5]/2
                Lv_matrix[1][0]=Lv[5]/2

                eps_matrix = self._get_eps_matrix(Lv.copy())

                #--- Calculating the cell_new matrix ---------------------------------------------------------
                i_matrix   = np.eye(3)
                def_matrix = i_matrix + eps_matrix
                cell_new      = np.dot(cell_old, def_matrix)

                if self.verbose:
                   box_center('raitio %.3f    Lv_matrix        eta_mtrix          def_matrix'%r)
                   prettyprint(np.hstack((Lv_matrix,eps_matrix,def_matrix)),precision=5)
                #------------------------------------------------------------------------------------------
                cont2 = cont2 + 1
                Defn_cont2  = os.path.join(Defn, Defn.split('/')[-1]+'_'+def_fmt2%cont2)
        
                def_structs_dict['Path'] = Defn_cont2
                def_structs_dict['eta']  = r
                def_structs_dict['Cell'] = cell_new.tolist()
                def_structs_list.append(def_structs_dict)

                #--- Writing the structure file -----------------------------------------------------------
                create_path(Defn_cont2,back=False)
                new_struct=Structure(cell_new,self.structure.species,self.structure.frac_coords)
                new_struct.to("POSCAR",os.path.join(Defn_cont2,Defn_cont2.split('/')[-1] +'.vasp'))
            def_params[def_key]={'Path':Defn,'LagrangianStrain':Ls_str[i],'DeformedCell':def_structs_list}
        #--------------------------------------------------------------------------------------------------
        dumpfn(def_params,os.path.join(self.workdir,'Deformed_Parameters.json'),indent=4)
        #os.chdir(pwd)
        self.to(os.path.join(self.workdir,felastic),option={"numb_points":numb_points,"max_lag_strain":max_lag_strain,"V0":V0})

#--------------------------------------------------------------------------------------------------
def post_elastic(args):
    '''
    post process function entry point
    '''
    logo()
    skip=args.skip
    dpi=args.dpi
    fmt=args.fmt
    order=args.order
    properties = args.properties
    workdir=os.path.join(os.getcwd(), args.properties+'_'+args.approach)
    elat=loadfn(os.path.join(workdir,felastic))
    if isinstance(elat,Elastic):
       pass
    else:
       elat=Elastic.from_dict(elat)
    # if not set, then the value of order will be reset
    if order==0:
       if elat.approach == 'stress':
          order = 3
       else:
          order = 4
    else:
       pass
    print(elat)
    if args.verbose:
       print('Default parameter for Elastic calculation post:')
       print(args)
       elat.verbose=args.verbose
       print('order ==> %d'%order)
    plot=args.plot

    if properties == 'elc':
       elat.calc_elastic_constant(poly_order=order,skip=skip,finput='input.yaml',plot=plot,fmt=fmt,dpi=dpi)
    else:
       elat.calc_stress_strain(skip=skip,finput='input.yaml',plot=plot,fmt=fmt,dpi=dpi)

#--------------------------------------------------------------------------------------------------
def init_elastic(args):
    '''
    init process function entry point
    '''
    logo()
    structure=args.config
    approach=args.approach
    numb_points =args.number
    properties = args.properties
    max_lag_strain =args.maxs
    direction=args.direction
    verbose=args.verbose
    back=args.back
    lag_strain_range=args.ranges
    workdir=os.path.join(os.getcwd(), args.properties+'_'+args.approach)
    elat=Elastic(structure,approach=approach,properties=properties,workdir=workdir,verbose=verbose)
    #if approach == 'stress':
    #   max_lag_strain = min([max_lag_strain,0.005])
    print(elat)
    if args.verbose:
       print('Default parameter for Elastic calculation init:')
       print(args)
    if properties == 'elc':
       elat.set_deformed_structure_for_elc(numb_points,max_lag_strain,back=back)
    else:
       if lag_strain_range:
          max_lag_strain=None
       elat.set_deformed_structure_for_ss(numb_points,max_lag_strain,lag_strain_range,direction,back=back)
        
#--------------------------------------------------------------------------------------------------
if __name__=="__main__":
   st=Structure.from_file('POSCAR')
   elat=Elastic(st,approach='energy',properties='elc',verbose=True)
   sts=elat.set_deformed_structure_for_elc(9,0.02,back=False) 
   #sts=elat.set_deformed_structure_for_elc(6,0.1,workdir='test1',back=True) 
   #sts=elat.set_deformed_structure_for_elc(6,0.4,back=False) 
   elat1=Elastic.from_dict(elat.as_dict())
   print(elat1)
   #print(elat1.calc_elastic_constant(skip=True))
   print(elat1.calc_elastic_constant(skip=False))
