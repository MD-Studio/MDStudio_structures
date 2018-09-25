# -*- coding: utf-8 -*-

"""
file: wamp_services.py

WAMP service methods the module exposes.
"""

from lie_structures.cheminfo_molhandle import (
     mol_addh, mol_attributes, mol_make3D, mol_read, mol_removeh, mol_write, mol_combine_rotations)


def create_path_file_obj(mol, extension='mol2'):
    """
    Encode the input files
    """
    return {
        'path': None, 'content': mol,
        'extension': extension}


class CheminfoMolhandleWampApi(object):
    """
    Cheminformatics molecule handling WAMP API
    """
    @staticmethod
    def read_mol(config):
        """Read molecular structure using `config` """

        return mol_read(
            config['mol']['content'], mol_format=config['mol']['extension'].lstrip('.'),
            toolkit=config['toolkit'])

    @staticmethod
    def get_output_format(config):
        """ Retrieve the format to store the output"""
        return config.get('output_format', config['mol']['extension'])

    def convert_structures(self, request, claims):
        """
        Convert input file format to a different format. For a detailed
        input description see the file:
           lie_structures/schemas/endpoints/convert_request_v1.json
        And for a detailed description of the output see:
           lie_structures/schemas/endpoints/convert_response_v1.json
        """
        molobject = self.read_mol(request)

        output = mol_write(molobject, mol_format=self.get_output_format(request), file_path=None)
        # Update session
        status = 'completed'

        return {'mol': create_path_file_obj(output), 'status': status}

    def addh_structures(self, request, claims):
        """
        Add hydrogens to the input structue. For a detailed
        input description see the file:
           lie_structures/schemas/endpoints/addh_request_v1.json
        And for a detailed description of the output see:
           lie_structures/schemaS/endpoints/addh_response_v1.json
        """

        molobject = mol_addh(
            self.read_mol(request),
            polaronly=request['polaronly'],
            correctForPH=request['correctForPH'],
            pH=request['pH'])

        output = mol_write(molobject, mol_format=self.get_output_format(request), file_path=None)

        # Update session
        status = 'completed'

        return {'mol': create_path_file_obj(output), 'status': status}

    def removeh_structures(self, request, claims):
        """
        Remove hydrogens from the input structure. For a detailed
        input description see the file:
           lie_structures/schemas/endpoints/removeh_request_v1.json
        And for a detailed description of the output see:
           lie_structures/schemas/endpoints/removeh_response_v1.json
        """
        molobject = mol_removeh(self.read_mol(request))

        output = mol_write(molobject, mol_format=self.get_output_format(request), file_path=None)

        # Update session
        status = 'completed'

        return {'mol': create_path_file_obj(output), 'status': status}

    def make3d_structures(self, request, claims):
        """
        Convert 1D or 2D structure representation to 3D.
        For a detailed
        input description see the file:
          lie_structures/schemas/endpoints/make3d_request_v1.json
        And for a detailed description of the output see:
          lie_structures/schemas/endpoints/make3d_response_v1.json
        """
        molobject = mol_make3D(
            self.read_mol(request),
            forcefield=request['forcefield'],
            localopt=request['localopt'],
            steps=request['steps'])

        output = mol_write(molobject, mol_format=self.get_output_format(request), file_path=None)

        # Update session
        status = 'completed'

        return {'mol': create_path_file_obj(output), 'status': status}

    def structure_attributes(self, request, claims):
        """
        Return common structure attributes
        For a detailed input description see the file:
          lie_structures/schemas/endpoints/info_request_v1.json

        And for a detailed description of the output see:
          lie_structures/schemas/endpoints/info_response_v1.json
        """
        # Retrieve the WAMP session information
        molobject = self.read_mol(request)
        attributes = mol_attributes(molobject) or {}

        # Update session
        status = 'completed'

        return {'status': status, 'attributes': attributes}

    def rotate_structures(self, request, claims):
        """
        Rotate the structure around an axis defined by x,y,z.
        For a detailed input description see the file:
          lie_structures/schemas/endpoints/rotate_request_v1.json

        And for a detailed description of the output see:
          lie_structures/schemas/endpoints/rotate_response_v1.json

        """
        # Read in the molecule
        molobject = self.read_mol(request)

        rotations = request['rotations']
        output = mol_combine_rotations(molobject, rotations=rotations)
        status = 'completed' if output is not None else 'failed'

        return {'status': status, 'mol': create_path_file_obj(output)}
