import os
import re
import numpy as np
from pymatgen.core import Structure

def vaspparser(finput="OUTCAR"):
    with open(finput) as rf:
        data = rf.read()
    nions = int(re.findall("NIONS\s*=\s*([0-9]*)", data)[0])
    natoms = np.array([int(x) for x in re.findall("ions per type\s*=\s*([\s0-9]*)", data)[0].split() ]   )
    lattice = np.array(
        [ x.split()[:3] for x in re.findall("direct lattice vectors.*\n([-+0-9.\s\n]*)length of vectors", data )[-1].split("\n")[:3]       ]
    ).astype(float)
    positions = np.array(
        [
            x.split()[:3]
            for x in re.findall(
                "position of ions in fractional coordinates.*\n([-+0-9.\s\n]*)",
                data,
            )[-1].split("\n")[: nions]
        ]
    ).astype(float)
    iontype = [int(x)  for x in re.findall("ions per type\s*=\s*([0-9\s]*)", data)[0].split()     ]
    _species = re.findall("VRHFIN\s*=([a-zA-Z\s]*):", data)
    if len(_species) ==0:
       print("Make sure that VRHFIN information is in the OUTCAR. Try to remove the NWRITE=1 in INCAR. \nYou can use 'grep VRHFIN POTCAR' to obtain the element information and add the output\nto your OUTCAR")
       os._exit(0)
    species=[]
    for n,el in zip(iontype,_species):
        species.extend([el]*n)
    structure = Structure(lattice, species, positions)
    # external pressure in kB -> GPa
    pressure = float(
        re.findall(r"external\s*pressure\s=\s*([-0-9.]*)\s*kB", data)[-1]
    )
    pressure = pressure / 10
    c = np.matrix(
        [
            x.split()[1:]
            for x in re.findall(
                "TOTAL ELASTIC MODULI \(kBar\).*\n.*\n.*\n([XYZ0-9.\s-]*)\n\s*-",
                data,
            )[0].split("\n")
        ][:6]
    ).astype(float)
    elastic_tensor = c.copy()
    for j in range(0, 6):
        elastic_tensor[3, j] = c[4, j]
        elastic_tensor[4, j] = c[5, j]
        elastic_tensor[5, j] = c[3, j]

    ctemp = elastic_tensor.copy()

    for i in range(0, 6):
        elastic_tensor[i, 3] = elastic_tensor[i, 4]
        elastic_tensor[i, 4] = elastic_tensor[i, 5]
        elastic_tensor[i, 5] = ctemp[i, 3]
    # Change the units of Cij from kBar to GPa
    elastic_tensor /= 10.0
    for i in range(0, 6):
        elastic_tensor[i, i] = elastic_tensor[i, i] - pressure
    elastic_tensor[0, 1] = elastic_tensor[0, 1] + pressure
    elastic_tensor[1, 0] = elastic_tensor[1, 0] + pressure
    elastic_tensor[0, 2] = elastic_tensor[0, 2] + pressure
    elastic_tensor[2, 0] = elastic_tensor[2, 0] + pressure
    elastic_tensor[1, 2] = elastic_tensor[1, 2] + pressure
    elastic_tensor[2, 1] = elastic_tensor[2, 1] + pressure
    C=elastic_tensor*structure.lattice.c/10 # GPa to N/m
    return structure, C.A

if __name__=='__main__':
   from utils import prettyprint
   st,C=vaspparser()
   print(st)
   print('elastic constant')
   prettyprint(C)
