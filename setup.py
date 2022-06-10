#!/usr/bin/env python
# -*- coding: utf-8 -*-

from os import path
from mech2d import NAME,SHORT_CMD
import setuptools, datetime

readme_file = path.join(path.dirname(path.abspath(__file__)), 'README.md')
try:
    from m2r import parse_from_file
    readme = parse_from_file(readme_file)
except ImportError:
    with open(readme_file) as f:
        readme = f.read()

today = datetime.date.today().strftime("%b-%d-%Y")
with open(path.join('mech2d', '_date.py'), 'w') as fp :
    fp.write('date = \'%s\'' % today)

install_requires=[
                  'pymatgen==2022.4.26',
                  'dpdispatcher==0.3.39',
                  'custodian==2021.2.8']

setuptools.setup(
    name=NAME,
    use_scm_version={'write_to': 'mech2d/_version.py'},
    setup_requires=['setuptools_scm'],
    author="Haidi Wang",
    author_email="haidi@hfut.edu.cn",
    description="mech2d: a calculator for two dimensional material mechanical properties",
    long_description=readme,
    long_description_content_type="text/markdown",
    url="https://gitee.com/haidi-hfut/mech2d",
    python_requires="~=3.6",
    packages=['mech2d', 
              'mech2d/calculation',
              'mech2d/scripts'
    ],
    data_files = [('mech2d/scripts/', ['mech2d/scripts/submit_slurm_serial.sh', ])],
    classifiers=[
        "Programming Language :: Python :: 3.6",
        "License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)",
    ],
    keywords='',
    install_requires=install_requires,    
        entry_points={
          'console_scripts': [
              SHORT_CMD+'= mech2d.main:main',
              'cvasp= mech2d.scripts.cvasp:main']
   }
)

