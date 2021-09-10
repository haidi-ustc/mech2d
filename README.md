# mech2d 

mech2d is a python package that used to calculate the mechanics properties of two dimenional materials, including elstaic constant tensor, stress-strain curve and other relavant properties. mech2d is user friendly to generate deformed structures, submit DFT calculation tasks and process results.

For more information, check the [document]()

## Installation 

You may clone the source code from gitee
```bash
git clone git@gitee.com:haidi-hfut/mech2d.git
cd mech2d 
python setup.py install
```
or directly download the source code zip package and then install it
```bash
wget https://gitee.com/haidi-hfut/mech2d/repository/archive/master.zip
unzip mech2d.zip
cd mech2d
python setup.py install
```

## Usage

The calculation of mechanical properties can be divided into three stages:

**init** , **run** and **post** 

all of operation in mech2d software is carried out by command line arguments, to get the help information by

```
m2d -h
```
it shows

```
usage: m2d [-h] [-v] {init,run,post} ...

Desctiption:
------------
mech2d is a convenient script that use to calculate the mechanical properties of
2D materials, including EOS, Stress-Strain Curve, elastic constants and revalant
properties. The script works based on several sub-commands with their own options.
To see the options for the sub-commands, type "m2d sub-command -h".

positional arguments:
  {init,run,post}
    init           Generating initial data for elastic systems.
    run            Run the DFT calculation for deformed structures.
    post           Post processing for elastic calculation.

optional arguments:
  -h, --help       show this help message and exit
  -v, --version    Display version
```

* the **init** stage used to generate the deformed structures. 

To obtain the deformed structures according to different symmetry, a relaxed structure if needed. Suppose you have a POSCAR file in current folder, using the following command to generate the deformed structures by stress fitting approach:
```bash
m2d init -c POSCAR -n  9 -m 0.02 -a stress
```
After running this command,  it will build a `elc_stess` folder with subfolders. 

```bash
elc_stress/
├── Def_1
│   ├── Def_1_001
│   ├── Def_1_002
│   ├── Def_1_003
│   ├── Def_1_004
│   ├── Def_1_005
│   ├── Def_1_006
│   ├── Def_1_007
│   ├── Def_1_008
│   └── Def_1_009
└── Def_2
    ├── Def_2_001
    ├── Def_2_002
    ├── Def_2_003
    ├── Def_2_004
    ├── Def_2_005
    ├── Def_2_006
    ├── Def_2_007
    ├── Def_2_008
    └── Def_2_009

```
likewise, the deformed structures based on  energy fitting approach can be obtained by:
```bash
m2d init -c POSCAR -n  9 -m 0.02 -a stress
```

the deformed structures used for stress strain curve calculation can be obtained by:
```bash
m2d init -c POSCAR -n 21 -r 0.0 0.2 -a stress -d 'xx' 'yy' -p ssc 
```
if the `-r` parameter exist, it means that Lagrangian stress will be set by a range, e.g. 0 0.2
if the `-m` parameter exist, it means the Lagrangian stress equals  [-max,+max]

> tips:  .vasp .xsf .cif format files are all supported

* the **run** stage used to conduct the DFT calculation. 

For this stage, the running command is :

```
m2d run -a stress input.yaml
```

the `input.yaml` use to set the  configure file, it includes machine , resources , tasks and code.

1. **machine** dict set the machine information, local or remote ? what kind of queue system will use ? what is the working directory
2. **resources** dict set the queue task infomation, how many cores will use? which queue to use?  and set the enviroment variables
2. **task** dict set the task related infomation, the excutable command , the input files (forward_files) and necessary result files that need to be copied back (back_files)
4. **code** dict set the DFT engine information. including code name, input file location 

An example for `input.yaml` is shown below:
```
---
machine:
  batch_type: Slurm
  context_type: LocalContext
  local_root: "./work"
  remote_root: "./work"
  remote_profile:
    hostname: localhost
    username: wang
    port: 22
    timeout: 10
resources:
  number_node: 1
  cpu_per_node: 48
  gpu_per_node: 0
  queue_name: batch
  task_max: 10
  group_size: 1
  custom_flags: 
   - "ulimit -s unlimited"
  module_list: 
   - "vasp/5.4.1"
  
  #source_list:
  # - "/opt/intel/parallel_studio_xe_2020.2.108/psxevars.sh intel64"
  #envs:
  #  PATH: "/opt/soft/vasp541:$PATH"
tasks:
  command: "mpirun -np 48 vasp_std"
  task_work_path:
  forward_files:
  - INCAR
  - KPOINTS
  - POTCAR
  - POSCAR
  backward_files:
  - runlog
  - errlog
  - OUTCAR
  - OSZICAR
  - vasprun.xml
  - CONTCAR
  outlog: runlog
  errlog: errlog
code:
  name: vasp
  input:
    INCAR: "./INCAR"
    #KPOINTS: "./KPOINTS"
    KPOINTS:
      kspacing: 5000
      kgamma: false
    POTCAR: "./POTCAR"
    #vdw_kernel: vdw_kernel.bindat
```

> The input.yaml can also be setted by JSON format, however it seems that
  Yaml format is much easier to read.


* the **post** stage used analysis the result and plot the result. 

For this stage, the running command is :

```
m2d post -a stress --plot
```

or for stress-strain calculation 

```
m2d post -a stress -p ssc --plot
```

