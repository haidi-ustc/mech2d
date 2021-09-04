#/usr/bin/env python
#-code- utf-8

import os
import sys
import time
import shutil
import numpy as np
from uuid import uuid4
from monty.json import  MSONable
from monty.serialization import loadfn,dumpfn
from pymatgen.core import Structure
from pymatgen.io.vasp import Incar,Kpoints,Poscar,Potcar, VaspInput, Vasprun
from mech2d.calculation.calculator import  Calculator

class VASP(Calculator,MSONable):
    name="VASP"
    def __init__(self,inputdir,structure,
                      incar=None,
                      kpoints=None,
                      potcar=None,
                      structure_fmt='POSCAR',
                      outputdir=None,
                      params=None
                   ):
        #print(params) 
        self.outputdir=outputdir if outputdir else inputdir 
        self.inputdir=inputdir
        self.structure_fmt=structure_fmt
        self._set_structure(structure)
        if params:
           incar=params['INCAR'] 
           kpoints=params['KPOINTS']
           potcar=params['POTCAR']
        self._set_incar(incar)
        self._set_kpoints(kpoints)
        self._set_potcar(potcar)
        if 'vdw_kernel' in params.keys():
            self._vdw=params['vdw_kernel']
        else:
            self._vdw=None

    def _set_structure(self,structure):
        if isinstance(structure,str):
           structure_file=os.path.join(self.inputdir,structure)
           if os.path.isfile(structure_file):
              self.structure=Structure.from_file(structure_file)
           else:
              self.structure=Structure.from_str(structure,fmt=self.structure_fmt)
        elif isinstance(structure,Structure):
           self.structure=structure
        elif isinstance(structure,dict):
           self.structure=Structure.from_dict(structure)
        else:
           raise RuntimeError("structure must be a file, Structure Str, dictor Structure obj.")

    def _set_incar(self,incar):

        if isinstance(incar,str):
           incar_file=os.path.join(self.inputdir,incar)
           if os.path.isfile(incar_file):
              self.incar=Incar.from_file(incar_file)
           elif os.path.isfile(incar):
              self.incar=Incar.from_file(incar)
           else:
              self.incar=Incar.from_string(incar)
        elif isinstance(incar,Incar):
           self.incar=incar
        elif isinstance(incar,dict):
           try:
               self.incar=Incar.from_dict(incar)
           except:
               self.incar=None
        elif incar is None:
             try:
                incar_file=os.path.join(self.inputdir,'INCAR')
                self.incar=Incar.from_file(incar_file)
             except:
                self.incar=None
        else:
           raise RuntimeError("incar must be a file, Incar Str, dictor Incar obj.")

    def _set_kpoints(self,kpoints):

        if isinstance(kpoints,str):
           kpoints_file=os.path.join(self.inputdir,kpoints)
           if os.path.isfile(kpoints_file):
              self.kpoints=Kpoints.from_file(kpoints_file)
           else:
              self.kpoints=Kpoints.from_file(kpoints)

        elif isinstance(kpoints,Kpoints):
           self.kpoints=kpoints

        elif isinstance(kpoints,dict):
           if 'kspacing' in kpoints.keys():
               if 'kgamma' in kpoints.keys():
                   self.kpoints=get_kpoints(self.structure, kpoints['kspacing'], kpoints['kgamma'], twod=True) 
               else:
                   self.kpoints=get_kpoints(self.structure, kpoints['kspacing'], False , twod=True) 
           else:
               self.kpoints=Kpoints.from_dict(kpoints)

        elif kpoints is None:
             try:
                kpoints_file=os.path.join(self.inputdir,"KPOINTS")
                self.kpoints=Kpoints.from_file(kpoints_file)
             except:
                self.kpoints=None
        else:
              raise RuntimeError("kpoints must be a file, Kpoints Str, dictor Kpoints obj.")

    def _set_potcar(self,potcar):
        if isinstance(potcar,str):
           potcar_file=os.path.join(self.inputdir,potcar)
           if os.path.isfile(potcar_file):
              self.potcar=Potcar.from_file(potcar_file)
           elif os.path.isfile(potcar):
              self.potcar=Potcar.from_file(potcar)
           else:
              raise RuntimeError("potcar must be a file, Potcar Str, dictor Potcar obj.")
        elif isinstance(potcar,Potcar):
           self.potcar=potcar
        elif isinstance(potcar,dict):
             if "functional" in potcar:
                 if  'symbol_set' in potcar:
                     self.potcar=Potcar(potcar['symbol_set'],functional=potcar["functional"])
                 else:
                     symbols=[el for el in self.structure.symbol_set]
                     self.potcar=Potcar(symbols,functional=potcar["functional"])
             else:    
                 try:
                     self.potcar=Potcar.from_dict(potcar)
                 except:
                      raise RuntimeError('unknow format for potcar dict')
        elif potcar is None:
             try:
                potcar_file=os.path.join(self.inputdir,'POTCAR')
                self.potcar=Potcar.from_file(potcar_file)
             except:
                raise RuntimeError("potcar must be a file, Potcar Str, dictor Potcar obj.")
        else:
             raise RuntimeError("potcar must be a file, Potcar Str, dictor Potcar obj.")

    def set_inputs(self):
        vi=VaspInput(self.incar,self.kpoints,Poscar(self.structure),self.potcar)
        vi.write_input(output_dir=self.outputdir)
        if self._vdw:
           shutil.copy(self._vdw,os.path.join(self.outputdir,'vdw_kernel.bindat'))

    @staticmethod
    def get_energy(outputdir):     
        try: 
           vr=Vasprun(os.path.join(outputdir,'vasprun.xml'))
           if vr.converged:
              return vr.final_energy
        except:
           return None

    @staticmethod
    def get_stress(outputdir):
        try: 
           vr=Vasprun(os.path.join(outputdir,'vasprun.xml'))
           if vr.converged:
              return vr.ionic_steps[-1]['stress']
        except:
           return None
    

def get_kpoints(struct, kspacing, kgamma, twod=True) :
    if isinstance(kspacing,list):
       assert len(kspacing)==3
    else:
       kspacing = [kspacing, kspacing, kspacing]

    rlat = struct.lattice.reciprocal_lattice.abc
    kpts = [max(1,int(np.ceil(2 * np.pi * ii / ks))) for ii,ks in zip(rlat,kspacing)]

    if twod:
       kpts[2]=1

    if kgamma: 
       return Kpoints.gamma_automatic(kpts=kpts)
    else:
       return Kpoints.monkhorst_automatic(kpts=kpts)
def test():
     print(VASP.get_energy('.') )
if __name__=="__main__":
     #test()
     d=loadfn('input.yaml')
     print(d)
     vasp=VASP('.', 'POSCAR', params= d['code']['input'])
     #kpoints = {'kspacing': 0.2,'kgamma':False} 
     #vasp=VASP('.','POSCAR',incar='INCAR',kpoints=kpoints,potcar='POTCAR') #{"functional":'PBE_54'})
     #vasp=VASP('.','POSCAR',incar='INCAR',kpoints=kpoints,potcar={"functional":'unvie_PBE_52'})
     print(vasp.name)
     print(vasp)
     print(vasp.incar)
     print(vasp.structure)
     #vasp.set_inputs()
