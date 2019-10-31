# MDStudio_structures

[![Build Status](https://travis-ci.org/MD-Studio/MDStudio_structures.svg?branch=master)](https://travis-ci.org/MD-Studio/MDStudio_structures)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/3c054785c5da46dfaad6dc3443d5653f)](https://www.codacy.com/manual/marcvdijk/MDStudio_structures?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=MD-Studio/MDStudio_structures&amp;utm_campaign=Badge_Grade)
[![codecov](https://codecov.io/gh/MD-Studio/MDStudio_structures/branch/master/graph/badge.svg)](https://codecov.io/gh/MD-Studio/MDStudio_structures)

![Configuration settings](mdstudio-logo.png)

Molecular structure based microservice for the MDStudio workflow environment.
This service combines the functionality of 9 popular cheminformatics packages through a common interface using 
[Cinfony](http://cinfony.github.io):

* webel
* silverwebel
* OpenBabel
* CDK
* RDkit
* Opsin
* Indigo
* JChem
* PyDPI

## Installation Quickstart
The best way to install and use MDStudio Structures and all of its dependencies is by running the precompiled
Docker container:
    
    docker pull mdstudio/mdstudio_structures