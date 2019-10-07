from __future__ import print_function

from mdstudio.deferred.chainable import chainable
from mdstudio.component.session import ComponentSession
from mdstudio.runner import main
from os.path import join

import numpy as np
import os
import pybel
import tempfile

file_path = os.path.realpath(__file__)
root = os.path.split(file_path)[0]
tempdir = tempfile.gettempdir()


def create_path_file_obj(path, extension='mol2'):
    """
    Encode the input files
    """
    return {
        u'path': path, u'content': read_file(path),
        u'extension': extension}


def create_smile_file_object(smile):
    """
    represent a molecule as a serialized file containing a smile
    """
    return {u'content': smile, u'extension': 'smi', u'path': None}


def read_file(fp):
    """
    Returns file content
    """
    with open(fp, 'r') as f:
        content = f.read()
    return content


def read_mol(mol, fmt='mol2'):
    """
    Read a molecular object either from a file or string
    """
    if os.path.isfile(mol):
        mol = pybel.readfile(fmt, mol).next()
    else:
        mol = pybel.readstring(fmt, mol)

    return mol


def compare_molecules(mol1, mol2):
    """
    Compare if the coordinates array are the same shape
    """
    m1 = read_mol(mol1)
    m2 = read_mol(mol2)

    arr = np.array([x.coords for x in m1])
    brr = np.array([x.coords for x in m2])

    return arr.shape == brr.shape


dict_convert = {
    u'output_format': u'mol2',
    u'workdir': tempdir,
    u'mol': create_smile_file_object(u'O1[C@@H](CCC1=O)CCC')
}

dict_make3d = {
    u'workdir': tempdir,
    u'output_format': u'mol2',
    u'mol': create_path_file_obj(join(root, 'files/structure.mol2'))
}

dict_addh = {
    u'workdir': tempdir,
    u'output_format': u'mol2',
    u'mol': create_path_file_obj(join(root, 'files/structure3D.mol2')),
    u'pH': 7.4,
    u'correctForPH': False
}

dict_info = {
    u'mol': create_path_file_obj(join(root, 'files/structure3D.mol2')),
    u'workdir': tempdir
}

dict_rotate = {
    u'workdir': tempdir,
    u'output_format': u'mol2',
    u'rotations': [
        [1, 0, 0, 90], [1, 0, 0, -90], [0, 1, 0, 90],
        [0, 1, 0, -90], [0, 0, 1, 90], [0, 0, 1, -90]],
    u'mol': create_path_file_obj(join(root, 'files/structureHs.mol2')),
}

dict_similarity = {
    u'mol_format': u'smi',
    u'ci_cutoff': 0.3617021276595745,
    u'workdir': tempdir,
    u'test_set': [create_smile_file_object(u'O1[C@@H](CCC1=O)CCC')],
    u'reference_set': map(create_smile_file_object, [
      u'c1(c(cccc1Nc1c(cccc1)C(=O)O)C)C',
      u'c12ccccc1nc1c(c2N)CCCC1',
      u'c1ccc(c(c1)[N+](=O)[O-])[C@H]1C(=C(NC(=C1C(=O)OC)C)C)C(=O)OC',
      u'c1cc(ccc1OCC)NC(=O)C',
      u'c12c3c(ccc1c(=O)cc(o2)c1ccccc1)cccc3',
      u'c1cc(cc(c1N/C=N/O)C)CCCC',
      u'c1(cccnc1Nc1cc(ccc1)C(F)(F)F)C(=O)O',
      u'c1ccc(c(c1C)OC[C@H](C)N)C',
      u'c1(OC[C@H](CNC(C)C)O)c2c(ccc1)cccc2',
      u'c12ccccc1cccc2',
      u'c12ccccc1cccc2C',
      u'c12ccccc1ccc(c2)C',
      u'c12ccccc1ccc(c2)F',
      u'c12ccccc1cc(cc2C)C',
      u'c12ccccc1c(ccc2C)C',
      u'c12cccc(c1cccc2Cl)Cl',
      u'c12cc(ccc1cc(cc2)C)C',
      u'C1CCC(=O)OCC1',
      u'O1[C@@H](CCC1=O)C',
      u'O1[C@@H](CCC1=O)CC',
      u'O1[C@@H](CCC1=O)CCC',
      u'O1[C@@H](CCC1=O)CCCCC',
      u'O1[C@@H](CCC1=O)CCCCCC',
      u'O1[C@@H](CCC1=O)CCCCCCC',
      u'C1C[C@@H](OC(=O)C1)CCCCC',
      u'c1c2c(ccc1)OC(=O)C2',
      u'c1c2c(ccc1)CC(=O)C2',
      u'c1c2c(ccc1)OCC2',
      u'c1c2c(ccc1)oc(=O)[nH]2',
      u'c1(ccccc1)c1ccccc1',
      u'c1c(cccc1)c1ccc(cc1)Cl',
      u'C1CCCC(C1)CCCC',
      u'[C@@H]1(OC(=O)CC1)c1ccccc1',
      u'c1(cc(oc(=O)c1)C)C',
      u'C1CC(=O)N([C@H]1c1cccnc1)C'
    ])
}


class RunStructures(ComponentSession):

    def authorize_request(self, uri, claims):
        return True

    @chainable
    def on_run(self):
        toolkits = yield self.call(
            'mdgroup.lie_structures.endpoint.supported_toolkits',
            {})
        if 'pybel' not in toolkits['toolkits']:
            raise AssertionError('PyBel toolkit not available')
        print('toolkits available: {}'.format(toolkits['toolkits']))

        convert = yield self.call(
            u'mdgroup.lie_structures.endpoint.convert', dict_convert)
        if not compare_molecules(convert['mol']['content'], join(root, u'files/structure.mol2')):
            raise AssertionError('Structures are not the same')
        print('converting {} from smile to mol2 succeeded!'.format(
            dict_convert['mol']))

        make3d = yield self.call(
            'mdgroup.lie_structures.endpoint.make3d', dict_make3d)
        if not compare_molecules(make3d['mol']['content'], join(root, 'files/structure3D.mol2')):
            raise AssertionError('Structures are not the same')
        print('successful creation of a 3D structure for {}'.format(
            dict_convert['mol']))

        addh = yield self.call(
            'mdgroup.lie_structures.endpoint.addh', dict_addh)
        if not compare_molecules(addh['mol']['content'], join(root, 'files/structureHs.mol2')):
            raise AssertionError('Structures are not the same')
        print('added hydrogens sucessfully!')

        info = yield self.call('mdgroup.lie_structures.endpoint.info', dict_info)
        print('info', info)
        atts = info['attributes']
        if not all((atts['formula'] == 'C7H12O2', atts['exactmass'] - 128.083729624 < 1e-5)):
            raise AssertionError('Structures are not the same')
        print('attributes information successfully retrieved!')

        rotate = yield self.call('mdgroup.lie_structures.endpoint.rotate', dict_rotate)
        if not compare_molecules(rotate['mol']['content'], join(root, 'files/rotations.mol2')):
            raise AssertionError('Structures are not the same')
        print('rotatation method succeeded!')

        similarity = yield self.call('mdgroup.lie_structures.endpoint.chemical_similarity', dict_similarity)
        if not similarity['results']['idx_max_sim']['0'] == 20:
            raise AssertionError('Similarity not the same')
        print('Similarity method succeeded!')


if __name__ == '__main__':
    main(RunStructures)
