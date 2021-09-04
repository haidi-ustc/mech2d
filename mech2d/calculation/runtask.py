#/usr/bin/env python
#-code- utf-8

import json
import os
from glob import glob
from pymatgen.core import Structure
from mech2d.calculation.vasp import VASP
from monty.serialization import loadfn,dumpfn
from dpdispatcher import Machine, Resources, Task, Submission

supported_calculators=["VASP"]

class RunTasks():
    def __init__(self, rootpath, structures, code):
        self.rootpath=rootpath
        self.structures = structures
        self.code = code
        calculator=code['name'].upper()
        inputdirs=[os.path.dirname(ii) for ii in structures] 
        task_dirs=['/'.join(ii.split('/')[-3:]) for ii in inputdirs]
        self.task_dirs=task_dirs

        assert calculator in supported_calculators
        tasks=[]
        if calculator=="VASP":
           for structure,inputdir in zip(structures,inputdirs):
               tasks.append(VASP(inputdir, structure, params=self.code['input']))
        elif calculator=="GPAW":
            raise RuntimeError("Unsupported calculators")
        else:
            raise RuntimeError("Unsupported calculators")
        self.tasks=tasks

    def preprocess(self):
        for task in self.tasks:
            task.set_inputs()

    def run(self,task_dict,machine,resources):
        task_list=[]
        for task in self.task_dirs:
            _task=Task(
                      command=task_dict['command'],
                      task_work_path=task,
                      forward_files=task_dict['forward_files'],
                      backward_files=task_dict['backward_files'],
                      outlog=task_dict['outlog'],
                      errlog=task_dict['errlog']
                   )
            task_list.append(_task)
        submission=Submission(work_base=self.rootpath,
                              machine=machine, 
                              resources=resources,
                              task_list=task_list,
                              forward_common_files=[], 
                              backward_common_files=[]
                            )
        submission.run_submission()


def run_elastic(args):
    assert os.path.exists(args.input)
    inputs=loadfn(args.input)
    keys=['machine','resources','tasks','code']
    for key in keys:
        assert key in inputs.keys()

    machine=Machine.load_from_dict(inputs['machine'])
    resources=Resources.load_from_dict(inputs['resources'])
    tasks=inputs['tasks']
    code=inputs['code']
    rootpath=os.getcwd()
    assert os.path.exists(os.path.join(rootpath,"Elastic2D.json"))
    with open('Elastic2D.json','r') as f:
         elatic2d=json.loads(f.read())
    structures=glob(os.path.join(elatic2d['outdir'],'Def_?','Def_?_???','Def_*vasp'))
    runtask=RunTasks(rootpath,structures,code)
    runtask.preprocess()
    runtask.run(tasks,machine,resources)

if __name__=='__main__':
   run_elastic()
