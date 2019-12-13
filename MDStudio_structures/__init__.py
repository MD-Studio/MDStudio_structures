# -*- coding: utf-8 -*-

"""
MDStudio Structures service for small molecule cheminformatics functions

Why import pandas here?:
Pandas and Openbabel share a common json library name
Adding the pandas import here avoids a REALLY NASTY race condition
If you remove the pandas import all the HELL may break loose
"""

import pandas
import os

from .cheminfo_pkgmanager import CinfonyPackageManager

__module__ = 'mdstudio_structures'
__docformat__ = 'restructuredtext'
__version__ = '{major:d}.{minor:d}'.format(major=1, minor=0)
__author__ = 'Marc van Dijk'
__status__ = 'release alpha1'
__date__ = '12 december 2019'
__licence__ = 'Apache Software License 2.0'
__url__ = 'https://github.com/MD-Studio/MDStudio_structures'
__copyright__ = "Copyright (c) VU University, Amsterdam"
__rootpath__ = os.path.dirname(__file__)
__all__ = ['toolkits']


# Load the toolkits
paths = {}
toolkits = CinfonyPackageManager(paths)
