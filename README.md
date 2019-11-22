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
MDStudio structures can be used in the MDStudio environment as Docker container or as standalone service.

### Install option 1 Pre-compiled Docker container
MDStudio structures can be installed quickly from a pre-compiled docker image hosted on DockerHub by:

    docker pull mdstudio/mdstudio_structures
    docker run (-d) mdstudio/mdstudio_structures

In this mode you will first need to launch the MDStudio environment itself in order for the MDStudio structures service
to connect to it. You can unify this behaviour by adding the MDStudio structures service to the MDStudio service 
environment as:

    MDStudio/docker-compose.yml:
        
        services:
           mdstudio_structures:
              image: mdstudio/mdstudio_structures
              links:
                - crossbar
              environment:
                - CROSSBAR_HOST=crossbar
              volumes:
                - ${WORKDIR}/mdstudio_structures:/tmp/mdstudio/mdstudio_structures

And optionally add `mdstudio_structures` to MDStudio/core/auth/settings.dev.yml for automatic authentication and 
authorization at startup.

### Install option 2 custom build Docker container
You can custom build the MDStudio structures Docker container by cloning the MDStudio_structures GitHub repository and run:

    docker build MDStudio_structures/ -t mdstudio/mdstudio_structures
    
After successful build of the container follow the steps starting from `docker run` in install option 1.

### Install option 3 standalone deployment of the service
If you prefer a custom installation over a (pre-)build docker container you can clone the MDStudio_structures GitHub
repository and install `mdstudio_structures` locally as:

    pip install (-e) mdstudio_structures/

Followed by:

    ./entry_point_mdstudio_structures.sh
    
or

    export MD_CONFIG_ENVIRONMENTS=dev,docker
    python -u -m mdstudio_structures