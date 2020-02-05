#!/usr/bin/env python
from setuptools import setup

setup(
    name='xenomapper2',
    version='2.0a1',
    author='Matthew Wakefield',
    author_email='matthew.wakefield@unimelb.edu.au',
    install_requires = [
      'setuptools >= 1.1.6',
      'pylazybam',
    ],
    packages=['xenomapper2',
              'xenomapper2.tests',
              ],
    zip_safe = True,
    include_package_data = True,
    url='https://git@bitbucket.org/genomematt/xenomapper2.git',
    license='BSD',
    entry_points={
        'console_scripts': ['xenomapper2 = pylazybam.cli:main',
                           ]
    },
    test_suite = "xenomapper2.tests.test_all",
    description='xenomapper2 - splitting mixed origin NGS reads',
    long_description='A tool for assessing reads mapped against two genomes'
                     'to identify uniquely mapping reads. Alt mapping aware'
                     'and suitable for uses such as MHC regions in GRCh38.',
    classifiers=[
          'Development Status :: 4 - Beta',
          'License :: OSI Approved :: BSD License',
          'Operating System :: POSIX',
          'Programming Language :: Python :: 3.6',
          'Programming Language :: Python :: 3.7',
          'Programming Language :: Python :: 3.8',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],

)
