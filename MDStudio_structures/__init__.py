# -*- coding: utf-8 -*-
"""
package:  MDStudio_structures
LIEStudio small molecule cheminformatics functions
"""
# Pandas and Openbabel share a common json library name
# Adding the pandas import here avoids a REALLY NASTY race condition
# If you remove the pandas import all the HELL may break loose
import pandas
import os

__licence__ = 'Apache Software License 2.0'
__rootpath__ = os.path.dirname(__file__)

from .cheminfo_pkgmanager import CinfonyPackageManager

# Load the toolkits
paths = {}
toolkits = CinfonyPackageManager(paths)

__all__ = ['toolkits']
