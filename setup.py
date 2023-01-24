#!/usr/bin/env python

from distutils.core import setup

setup(
    name='RREFinder',
    scripts=['RRE.py', 'download_RRE_databases.py', 'setup_RRE_exploratory.py'],
    description='Bioinformatic application for the detection of RREs in protein sequences of interest',
    url='https://github.com/Alexamk/RREFinder',
    author='Alexander Kloosterman',
    author_email='alexander.kloosterman@ki.se',
    )
