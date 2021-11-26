#!/usr/bin/env python
# coding: utf-8
# Copyright (c) 2021.

import argparse
import sys
import itertools
from mech2d.mechanics import init_elastic,post_elastic
from mech2d.calculation.runtask import run_elastic
from mech2d import NAME

__author__ = "Haidi Wang"
__copyright__ = "Copyright 2021"
__maintainer__ = "Haidi Wang"
__email__ = "haidi@hfut.edu.cn"

def get_version():
  try:
     from mech2d._version import version
  except:
     version="Unknow" 
  return version 

def main():
    parser = argparse.ArgumentParser(prog='m2d',
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description="""Desctiption:\n------------
mech2d is a convenient script that use to calculate the mechanical properties of
2D materials, including Stress-Strain Curve, elastic constants and revalant
properties. The script works based on several sub-commands with their own options.
To see the options for the sub-commands, type "m2d sub-command -h".""")
    parser.add_argument('-v', '--version', action='version', version=get_version(),help='Display version')

    subparsers = parser.add_subparsers()

    #-------------
    # init 
    parser_init = subparsers.add_parser(
        "init", help="Generating initial data for elastic systems.")
    parser_init.add_argument('-c','--config', type=str, default='POSCAR', help="The structure filename. Supported format: ['.vasp','POSCAR','.cif','.xsf']")
    parser_init.add_argument('-a','--approach', type=str, default='energy', choices=['stress', 'energy'],help="Support 'Energy' or 'Stress' approach.")
    parser_init.add_argument('-m','--maxs', type=float, default=0.05, help="For elastic constant calculation, it stands for the maximum Lagrangian strain, suggested value is [0.030, 0.150] for Energy approach, [0.0010, 0.0050] for Stress approach; for stress strain cuver calcuation, this value has no above limitation")
    parser_init.add_argument('-n','--number', type=int, default=5, help="The number of the deformed structures [odd number > 4].")
    #parser_init.add_argument('-o','--outdir', type=str, default=None, help="The output directory for the deformed structures.")
    parser_init.add_argument('-d','--direction', type=str, default=['xx'],nargs='+', choices=['xx','yy','bi','xy'], help="The direction used for stress strain curve,  default value: 'xx'. 'xx' for 'x' direction; 'yy' for 'y' direction; 'bi' for bi-Axis strain and 'xy' for shear strain.")
    parser_init.add_argument('-r','--ranges', type=float, default=None, nargs='+', help="The Lagrangian strain range used for stress-strain curve calculation. e.g. 0.0 0.2 ")
    parser_init.add_argument('-p','--properties', type=str, default='elc', choices=['elc', 'ssc'], help="What do you want to calcuation? elastic constant or stress strain curve? default value: 'elc'.")
    parser_init.add_argument('-v','--verbose', action="store_true", help="print verbose information or not.")
    parser_init.add_argument('-b','--back', action="store_true" , help="Whether back the old folder? default value: False.")

    parser_init.set_defaults(func=init_elastic)
    
    #-------------
    #run 
    parser_run= subparsers.add_parser(
        "run", help="Run the DFT calculation for deformed structures.")
    parser_run.add_argument('-a','--approach', type=str, default='energy', choices=['stress', 'energy'],help="Support 'Energy' or 'Stress' approach.")
    parser_run.add_argument('-p','--properties', type=str, default='elc', choices=['elc', 'ssc'], help="What do you want to calcuation? elastic constant or stress strain curve? default value: 'elc'.")
    parser_run.add_argument("input", type=str,help="input file for supplying information about DFT calculation, json/yaml format. The 'machine', 'tasks', 'code', 'resources' should be supplied.")
    parser_run.add_argument('--manual', action="store_true", help="manual model, only for generating the input files without runing")
    parser_run.add_argument('-v','--verbose', action="store_true", help="print verbose information or not.")
    parser_run.set_defaults(func=run_elastic)
    
    #-------------
    #post
    parser_post = subparsers.add_parser(
        "post", help="Post processing for elastic calculation.")
    parser_post.add_argument('-a','--approach', type=str, default='energy', choices=['stress', 'energy'],help="Support 'Energy' or 'Stress' approach.")
    parser_post.add_argument('-p','--properties', type=str, default='elc', choices=['elc', 'ssc'], help="What do you want to calcuation? elastic constant or stress strain curve? default value: 'elc'.")
    parser_post.add_argument('--skip', action="store_true" , help="Whether skip the data parsing ? if true, it means the Def_*_Energy.dat should be exists in corresponding folder. default value: False.")
    parser_post.add_argument('-o','--order', type=int, default=0, help="The order of polynomial for fitting. Default value: 4 for strain-stress approach and 3 for stress-strain method")
    parser_post.add_argument('-f','--fmt', type=str, default='png', help="The format of output figure. Default value: .jpg")
    parser_post.add_argument('-d','--dpi', type=int, default=100, help="The resolution of output figure. Default value: 100")
    parser_post.add_argument('--plot', action="store_true", help="plot the figures")
    parser_post.add_argument('-v','--verbose', action="store_true", help="print verbose information or not.")
    parser_post.set_defaults(func=post_elastic)

    try:
        import argcomplete
        argcomplete.autocomplete(parser)
    except ImportError:
        # argcomplete not present.
        pass

    args = parser.parse_args()

    try:
        getattr(args, "func")
    except AttributeError:
        parser.print_help()
        sys.exit(0)
    args.func(args)

if __name__ == "__main__":
    main() 
