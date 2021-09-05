#!/usr/bin/env python
# coding: utf-8
# Copyright (c) 2021.

import argparse
import sys
import itertools
from mech2d.mechanics import init_elastic,post_elastic
from mech2d.calculation.runtask import run_elastic
from mech2d.plot import plot_elastic
from mech2d import NAME

def get_version():
  try:
     from mech2d._version import version
  except:
     version="Unknow" 
  return version 

__author__ = "Haidi Wang"
__copyright__ = "Copyright 2021"
__maintainer__ = "Haidi Wang"
__email__ = ""


def main():
    parser = argparse.ArgumentParser(prog='m2d',
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description="""Desctiption:\n------------
mech2d is a convenient script that use to calculate the mechanical properties of
2D materials, including EOS, Stress-Strain Curve, elastic constants and revalant
properties. The script works based on several sub-commands with their own options.
To see the options for the sub-commands, type "m2d sub-command -h".""")
    parser.add_argument('-v', '--version', action='version', version=get_version(),help='Display version')

    subparsers = parser.add_subparsers()

    #-------------
    # init 
    parser_init = subparsers.add_parser(
        "init", help="Generating initial data for elastic systems.")
    parser_init.add_argument('-c','--config', type=str, default='POSCAR', help="The structure filename. Supported format: ['.vasp','POSCAR','.cif','.xsf']")
    parser_init.add_argument('-s','--strategy', type=str, default='energy', choices=['stress', 'energy'],help="Support 'Energy' or 'Stress' strategy.")
    parser_init.add_argument('-m','--maxs', type=float, default=0.05, help="For elastic constant calculation, it stands for the maximum Lagrangian strain, suggested value is [0.030, 0.150] for Energy strategy, [0.0010, 0.0050] for Stress strategy; for stress strain cuver calcuation, this value has no above limitation")
    parser_init.add_argument('-n','--number', type=int, default=5, help="The number of the deformed structures [odd number > 4].")
    #parser_init.add_argument('-o','--outdir', type=str, default=None, help="The output directory for the deformed structures.")
    parser_init.add_argument('-d','--direction', type=str, default='x', choices=['xx','yy','bi','xy'], help="The direction used for stress strain curve,  default value: 'xx'. 'xx' for 'x' direction; 'yy' for 'y' direction; 'bi' for bi-axis strain and 'xy' for shear strain.")
    parser_init.add_argument('-p','--properties', type=str, default='elc', choices=['elc', 'ssc'], help="What do you want to calcuation? elastic constant or stress strain curve? default value: 'elc'.")
    parser_init.add_argument('-v','--verbose', action="store_true", help="print verbose information or not.")
    parser_init.add_argument('-b','--back', action="store_true" , help="Whether back the old folder? default value: False.")

    parser_init.set_defaults(func=init_elastic)
    
    #-------------
    #run 
    parser_run= subparsers.add_parser(
        "run", help="Run the DFT calculation for deformed structures.")
    parser_run.add_argument('-s','--strategy', type=str, default='energy', choices=['stress', 'energy'],help="Support 'Energy' or 'Stress' strategy.")
    parser_run.add_argument('-p','--properties', type=str, default='elc', choices=['elc', 'ssc'], help="What do you want to calcuation? elastic constant or stress strain curve? default value: 'elc'.")
    parser_run.add_argument("input", type=str,help="input file for supplying information about DFT calculation, json/yaml format. The 'machine', 'tasks', 'code', 'resources' should be supplied.")
    parser_run.add_argument('-v','--verbose', action="store_true", help="print verbose information or not.")
    parser_run.set_defaults(func=run_elastic)
    
    #-------------
    #post
    parser_post = subparsers.add_parser(
        "post", help="Post processing for elastic calculation.")
    parser_post.add_argument('-s','--strategy', type=str, default='energy', choices=['stress', 'energy'],help="Support 'Energy' or 'Stress' strategy.")
    parser_post.add_argument('-p','--properties', type=str, default='elc', choices=['elc', 'ssc'], help="What do you want to calcuation? elastic constant or stress strain curve? default value: 'elc'.")
    parser_post.add_argument('--skip', action="store_true" , help="Whether skip the data parsing ? if true, it means the Def_*_Energy.dat should be exists in corresponding folder. default value: False.")
    parser_post.add_argument('-o','--order', type=int, default=4, help="The order of polynomial for fitting. Default value: 4")
    parser_post.add_argument('-v','--verbose', action="store_true", help="print verbose information or not.")
    parser_post.set_defaults(func=post_elastic)

    #-------------
    #plot
    parser_plot = subparsers.add_parser(
        "plot", help="Plot the results")
    parser_plot.add_argument('-f','--filename', type=str, default='mech2d', help="The filename of output figure.")
    parser_plot.add_argument('-d','--datafile', type=str, default='result.json', help="The data file for plotting.")
    parser_plot.add_argument('-p','--properties', type=str, default='elc', choices=['elc', 'ssc'], help="What do you want to plot? elastic constant or stress strain curve? default value: 'elc'.")
    parser_plot.add_argument('-v','--verbose', action="store_true", help="print verbose information or not.")
    parser_plot.set_defaults(func=plot_elastic)

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
