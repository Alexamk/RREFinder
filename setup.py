#!/usr/bin/env python

import os
from distutils.core import setup

def generate_data_files():
    data_files = []
    for directory in ['data','lib','sample_input']:
        for path, dirs, files in os.walk(directory):
            install_dir = os.path.join('bin',path)
            list_entry = (install_dir, [os.path.join(path, f) for f in files ])
            data_files.append(list_entry)
    data_files.append(('bin',['config.ini']))
    return data_files

setup(
    name='RREFinder',
    version='v1.0.1',
    py_modules=['RRE','download_RRE_databases','setup_RRE_exploratory'],
    scripts=['RRE.py','download_RRE_databases.py','setup_RRE_exploratory.py'],
    data_files=generate_data_files(),
    description='Bioinformatic application for the detection of RREs in protein sequences of interest',
    url='https://github.com/Alexamk/RREFinder',
    author='Alexander Kloosterman',
    author_email='alexander.kloosterman@ki.se',
    license='GNU Affero v3'
    )
