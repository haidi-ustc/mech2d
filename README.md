# mech2d 

mech2d is a python package that used to calculate the mechanics properties of two dimenional materials, including elstaic constant tensor, stress-strain curve and other relavant properties. mech2d is user friendly to generate deformed structures, submit DFT calculation tasks and process results.

For more information, check the [document](https://mech2d.readthedocs.io/en/latest/) or [paper](https://www.preprints.org/manuscript/202304.0676/v1)

## Installation 

You may clone the source code from gitee
```bash
git clone git@gitee.com:haidi-hfut/mech2d.git
```

Before you install the mech2d. It is better to create a virtual python enviroment via conda

```bash
conda create -n mech2d python=3.10
```

After creating the virtual enviroment,  activate it and install the mech2d

```
conda activate mech2d
cd mech2d 
pip install .
```
or directly download the source code zip package and then install it
```bash
wget https://gitee.com/haidi-hfut/mech2d/repository/archive/master.zip
unzip mech2d.zip
cd mech2d
pip install .
```
## Method

Please refer to [mech2d](https://doi.org/10.3390/molecules28114337)
The polar plot of Young's modulus and Poisson's ratio is obtained by following equation:
```
v_{zz} & = \frac{C_{12}}{C_{22}} \\
d_1 & = \frac{C_{11}}{C_{22}} + 1 - \frac{C_{11} C_{22} - C_{12}^2}{C_{22} C_{66}} \\
d_2 & = -\left(2 \frac{C_{12}}{C_{22}} - \frac{C_{11} C_{22} - C_{12}^2}{C_{22} C_{66}}\right) \\
d_3 & = \frac{C_{11}}{C_{22}} \\
Y_{zz} & = \frac{C_{11} C_{22} - C_{12}^2}{C_{22}} \\
\theta & \in [0, 2\pi] \text{ with 360 points} \\
E(\theta) & = \frac{Y_{zz}}{\cos(\theta)^4 + d_2 \cos(\theta)^2 \sin(\theta)^2 + d_3 \sin(\theta)^4} \\
V(\theta) & = \frac{v_{zz} \cos(\theta)^4 - d_1 \cos(\theta)^2 \sin(\theta)^2 + v_{zz} \sin(\theta)^4}{\cos(\theta)^4 + d_2 \cos(\theta)^2 \sin(\theta)^2 + d_3 \sin(\theta)^4}
\en
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

use subcomand to show the help information
```
m2d init -h
```
it shows 
```
usage: m2d init [-h] [-c CONFIG] [-a {stress,energy}] [-m MAXS] [-n NUMBER] [-d {xx,yy,bi,xy} [{xx,yy,bi,xy} ...]] [-r RANGES [RANGES ...]]
                [-p {elc,ssc}] [-v] [-b]

optional arguments:
  -h, --help            show this help message and exit
  -c CONFIG, --config CONFIG
                        The structure filename. Supported format: ['.vasp','POSCAR','.cif','.xsf']
  -a {stress,energy}, --approach {stress,energy}
                        Support 'Energy' or 'Stress' approach.
  -m MAXS, --maxs MAXS  For elastic constant calculation, it stands for the maximum Lagrangian strain, suggested value is [0.030, 0.150] for
                        Energy approach, [0.0010, 0.0050] for Stress approach; for stress strain cuver calcuation, this value has no above
                        limitation
  -n NUMBER, --number NUMBER
                        The number of the deformed structures [odd number > 4].
  -d {xx,yy,bi,xy} [{xx,yy,bi,xy} ...], --direction {xx,yy,bi,xy} [{xx,yy,bi,xy} ...]
                        The direction used for stress strain curve, default value: 'xx'. 'xx' for 'x' direction; 'yy' for 'y' direction; 'bi'
                        for bi-axis strain and 'xy' for shear strain.
  -r RANGES [RANGES ...], --ranges RANGES [RANGES ...]
                        The Lagrangian strain range used for stress-strain curve calculation. e.g. 0.0 0.2
  -p {elc,ssc}, --properties {elc,ssc}
                        What do you want to calcuation? elastic constant or stress strain curve? default value: 'elc'.
  -v, --verbose         print verbose information or not.
  -b, --back            Whether back the old folder? default value: False.
```

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

use subcomand to show the help information
```
m2d run -h
```
it shows 
```
usage: m2d run [-h] [-a {stress,energy}] [-p {elc,ssc}] [--manual] [-v] input

positional arguments:
  input                 input file for supplying information about DFT
                        calculation, json/yaml format. The 'machine', 'tasks',
                        'code', 'resources' should be supplied.

optional arguments:
  -h, --help            show this help message and exit
  -a {stress,energy}, --approach {stress,energy}
                        Support 'Energy' or 'Stress' approach.
  -p {elc,ssc}, --properties {elc,ssc}
                        What do you want to calcuation? elastic constant or
                        stress strain curve? default value: 'elc'.
  --manual              manual model, only for generating the input files
                        without runing
  -v, --verbose         print verbose information or not.
```

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

use subcomand to show the help information
```
m2d post -h
```
it shows 
```
usage: m2d post [-h] [-a {stress,energy}] [-i INPUTFILE] [-p {elc,ssc}]
                [--skip] [-o ORDER] [-f FMT] [-d DPI] [--plot] [-v]

optional arguments:
  -h, --help            show this help message and exit
  -a {stress,energy}, --approach {stress,energy}
                        Support 'Energy' or 'Stress' approach.
  -i INPUTFILE, --inputfile INPUTFILE
                        Parsing elastic constant tensor from input file
  -p {elc,ssc}, --properties {elc,ssc}
                        What do you want to calcuation? elastic constant or
                        stress strain curve? default value: 'elc'.
  --skip                Whether skip the data parsing ? if true, it means the
                        Def_*_Energy.dat should be exists in corresponding
                        folder. default value: False.
  -o ORDER, --order ORDER
                        The order of polynomial for fitting. Default value: 4
                        for strain-stress approach and 3 for stress-strain
                        method
  -f FMT, --fmt FMT     The format of output figure. Default value: .jpg
  -d DPI, --dpi DPI     The resolution of output figure. Default value: 100
  --plot                plot the figures
  -v, --verbose         print verbose information or not.
```

For this stage, the running command is :

```
m2d post -a stress --plot
```
A figure named stress-EV.jpg will be generated.

![EV](https://note.youdao.com/yws/api/personal/file/A22D69ECC94C45138EBBB385DC4AAB95?method=download&shareKey=eebde950bb2defb3f61c1ec285a50b2f)

or for stress-strain calculation 

```
m2d post -a stress -p ssc --plot
```

## How to Cite

Wang, H.; Li, T.; Liu, X.; Zhu, W.; Chen, Z.; Li, Z.; Yang, J. mech2d: An Efficient Tool for High-Throughput Calculation of Mechanical Properties for Two-Dimensional Materials. [Link](https://doi.org/10.3390/molecules28114337)
