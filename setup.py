#! /usr/bin/env python
# -*- coding: utf-8 -*-

# package: lie_docking
# file: setup.py
#
# Part of ‘lie_docking’, a package providing molecular docking functionality
# for the LIEStudio package.
#
# Copyright © 2016 Marc van Dijk, VU University Amsterdam, the Netherlands
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from setuptools import setup, find_packages

distribution_name = 'lie_structures'

setup(
    name=distribution_name,
    version=0.2,
    description='MDStudio structure database module',
    author="""
    Marc van Dijk - VU University - Amsterdam
    Paul Visscher - Zefiros Software (www.zefiros.eu)
    Felipe Zapata - eScience Center (https://www.esciencecenter.nl/)""",
    author_email=['m4.van.dijk@vu.nl', 'f.zapata@esciencecenter.nl'],
    url='https://github.com/MD-Studio/MDStudio',
    license='Apache Software License 2.0',
    keywords='MDStudio structures database',
    platforms=['Any'],
    packages=find_packages(),
    py_modules=[distribution_name],
    install_requires=['biopython', 'cinfony==1.2', 'openbabel', 'pandas',
                      'pydpi', 'JPype1', 'Pillow', 'retrying', 'scipy'],
    extras_require={'test': ['unittest2']},
    dependency_links=["https://github.com/cinfony/cinfony/tarball/master#egg=cinfony-1.2"],
    include_package_data=True,
    zip_safe=True,
    classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: Apache Software License',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Operating System :: OS Independent',
        'Intended Audience :: Science/Research',
    ],
)
