#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import sys
import argparse
import numpy as np
from pymatgen.io.vasp import Kpoints
from custodian.vasp.jobs import VaspJob
from custodian.vasp.validators import VaspFilesValidator,VasprunXMLValidator
from custodian.vasp.handlers import VaspErrorHandler,UnconvergedErrorHandler, \
        NonConvergingErrorHandler,FrozenJobErrorHandler,StdErrHandler,\
        WalltimeHandler,PositiveEnergyErrorHandler
from custodian import Custodian

handlers=[VaspErrorHandler(),FrozenJobErrorHandler(),StdErrHandler(),NonConvergingErrorHandler(),
           WalltimeHandler(),PositiveEnergyErrorHandler(),UnconvergedErrorHandler()]
validators=[VaspFilesValidator(),VasprunXMLValidator()]

def relaxation_static_run(vasp_cmd,half_kpts_first_relax=False,auto_npar=True,backup=False,auto_continue=False):
    incar_update = {"ISTART": 1,"NSW":0}
    settings_overide_1 = None
    settings_overide_2 = [
        {"dict": "INCAR", "action": {"_set": incar_update}},
        {"file": "CONTCAR", "action": {"_file_copy": {"dest": "POSCAR"}}},
    ]
    if half_kpts_first_relax and os.path.exists("KPOINTS") and os.path.exists("POSCAR"):
        kpts = Kpoints.from_file("KPOINTS")
        orig_kpts_dict = kpts.as_dict()
        # lattice vectors with length < 8 will get >1 KPOINT
        kpts.kpts = np.round(np.maximum(np.array(kpts.kpts) / 2, 1)).astype(int).tolist()
        low_kpts_dict = kpts.as_dict()
        settings_overide_1 = [{"dict": "KPOINTS", "action": {"_set": low_kpts_dict}}]
        settings_overide_2.append({"dict": "KPOINTS", "action": {"_set": orig_kpts_dict}})

        return [
            VaspJob(
                vasp_cmd,
                final=False,
                suffix=".relax",
                backup=backup,
                auto_npar=auto_npar,
                auto_continue=auto_continue,
                settings_override=settings_overide_1,
            ),
            VaspJob(
                vasp_cmd,
                final=True,
                backup=backup,
                suffix="",
                auto_npar=auto_npar,
                auto_continue=auto_continue,
                settings_override=settings_overide_2,
            ),
        ]

def runvasp(cmd,max_errors=3,half_kpts_first_relax=False,auto_npar=True):
  
    """
    cmd example:
    cmd=['mpirun', '-np', '48' , '/opt/soft/vasp541/vasp_fix_z_intel2015']
    """
    jobs=relaxation_static_run(cmd,half_kpts_first_relax=half_kpts_first_relax,auto_npar=auto_npar)
    c = Custodian(handlers, jobs, validators=validators,max_errors=max_errors)
    c.run()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("cmd", type=str,
                        help="""The command for runing vasp, e.g.,
                              'mpirun -np 32 /path/vasp_std' or
                              'srun /path/vasp_std'
                             """)
    parser.add_argument("maxerrors", type=int,
                        help="The maximum error time for runing vasp")

    args = parser.parse_args()
    cmd=args.cmd.split()
    max_errors=args.maxerrors
    runvasp(cmd=cmd,max_errors=max_errors)

if __name__=='__main__':
  main()
  # vasp="/sharedext4/vasp/vasp.5.4.4/bin/vasp_std"
  # ncpu="48"
  # fcmd="/opt/soft/vasp541/vasp_fix_z_intel2015"
  # max_errors=3
  # runvasp(cmd=['mpirun', '-np', ncpu, fcmd],max_errors=max_errors)
