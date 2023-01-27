#!/usr/bin/env python

import os
from distutils.core import setup

def generate_package_data():
    data_files = []
    for directory in ['data','lib','sample_input']:
        for path, dirs, files in os.walk(directory):
            for f in files:
                data_files.append(os.path.join(path, f))
    data_files.append('config.ini')
    package_data={'RREFinder': data_files}
    return package_data

setup(
    name='RREFinder',
    version='v1.0.1',
    py_modules=['RRE','download_RRE_databases','setup_RRE_exploratory'],
    scripts=['RRE.py','download_RRE_databases.py','setup_RRE_exploratory.py'],
    packages=['RREFinder'],
    package_dir={'RREFinder': '.'},
    package_data=generate_package_data(),
    description='Bioinformatic application for the detection of RREs in protein sequences of interest',
    url='https://github.com/Alexamk/RREFinder',
    author='Alexander Kloosterman',
    author_email='alexander.kloosterman@ki.se',
    license='GNU Affero v3'
    )
