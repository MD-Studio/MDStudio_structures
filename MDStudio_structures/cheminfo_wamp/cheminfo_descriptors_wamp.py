# -*- coding: utf-8 -*-

from MDStudio_structures.cheminfo_molhandle import mol_read, mol_validate_file_object


class CheminfoDescriptorsWampApi(object):
    """
    Cheminformatics descriptors WAMP API
    """

    @staticmethod
    def get_descriptors(request, claims):

        # Import the molecule
        mol = mol_validate_file_object(request['mol'])
        molobject = mol_read(mol['content'], mol_format=mol['extension'], toolkit=request["toolkit"])
        desc = molobject.calcdesc()

        if desc is not None:
            status = 'completed'
            output = desc
        else:
            status = 'failed'
            output = None

        return {'session': status, 'descriptors': output}
