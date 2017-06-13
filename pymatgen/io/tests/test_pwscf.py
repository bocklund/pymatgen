# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

__author__ = 'Shyue Ping Ong'
__copyright__ = 'Copyright 2013, The Materials Project'
__version__ = '0.1'
__maintainer__ = 'Shyue Ping Ong'
__email__ = 'ongsp@ucsd.edu'
__date__ = '3/28/15'

import unittest
import os

from pymatgen.io.pwscf import PWInput, PWInputError, PWOutput
from pymatgen.util.testing import PymatgenTest

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        'test_files')

class PWInputTest(PymatgenTest):

    def test_init(self):
        s = self.get_structure("Li2O")
        self.assertRaises(
            PWInputError, PWInput,
            s,
            control={"calculation": "scf", "pseudo_dir": './'},
            pseudo={"Li": "Li.pbe-n-kjpaw_psl.0.1.UPF"}
        )

    def test_str(self):
        s = self.get_structure("Li2O")

        pw = PWInput(s,
                     control={"calculation": "scf", "pseudo_dir": './'},
                     pseudo={"Li": "Li.pbe-n-kjpaw_psl.0.1.UPF",
                             "O": "O.pbe-n-kjpaw_psl.0.1.UPF"},
                     system={"ecutwfc": 50})
        ans = """&CONTROL
  calculation = 'scf',
  pseudo_dir = './',
/
&SYSTEM
  ecutwfc = 50,
  ibrav = 0,
  nat = 3,
  ntyp = 2,
/
&ELECTRONS
/
&IONS
/
&CELL
/
ATOMIC_SPECIES
  Li  6.9410 Li.pbe-n-kjpaw_psl.0.1.UPF
  O  15.9994 O.pbe-n-kjpaw_psl.0.1.UPF
ATOMIC_POSITIONS crystal
  Li 0.250000 0.250000 0.250000
  Li 0.750000 0.750000 0.750000
  O 0.000000 0.000000 0.000000
K_POINTS automatic
  1 1 1 0 0 0
CELL_PARAMETERS angstrom
  -2.305000 -2.305000 0.000000
  -2.305000 0.000000 -2.305000
  0.000000 -2.305000 -2.305000
"""
        self.assertEqual(pw.__str__().strip(), ans.strip())

    def test_read_write_equivalence(self):
        """PWInput can read the structure it writes
        """
        s = self.get_structure("Li2O")

        pw = PWInput(s,
                     control={"calculation": "scf", "pseudo_dir": './'},
                     pseudo={"Li": "Li.pbe-n-kjpaw_psl.0.1.UPF",
                             "O": "O.pbe-n-kjpaw_psl.0.1.UPF"},
                     system={"ecutwfc": 50})
        ans = """&CONTROL
  calculation = 'scf',
  pseudo_dir = './',
/
&SYSTEM
  ecutwfc = 50,
  ibrav = 0,
  nat = 3,
  ntyp = 2,
/
&ELECTRONS
/
&IONS
/
&CELL
/
ATOMIC_SPECIES
  Li  6.9410 Li.pbe-n-kjpaw_psl.0.1.UPF
  O  15.9994 O.pbe-n-kjpaw_psl.0.1.UPF
ATOMIC_POSITIONS crystal
  Li 0.250000 0.250000 0.250000
  Li 0.750000 0.750000 0.750000
  O 0.000000 0.000000 0.000000
K_POINTS automatic
  1 1 1 0 0 0
CELL_PARAMETERS angstrom
  -2.305000 -2.305000 0.000000
  -2.305000 0.000000 -2.305000
  0.000000 -2.305000 -2.305000
"""
        pw_read = PWInput.from_string(ans)
        # we have to have a small custom version of __eq__
        self.assertIsNotNone(pw_read.structure, msg='Structure read from string is None')
        self.assertEqual(len(pw.structure), len(pw_read.structure),  msg='Structures do not have the same number of sites')
        # Handle float precision and original and read structures
        self.assertArrayAlmostEqual(pw.structure.lattice.matrix,pw_read.structure.lattice.matrix, err_msg='Structures are not equal within 7 decimal places')
        # We can't get site properties like charge from the input file, so we just check that the elements are there
        for specie in pw.structure.species:
            self.assertIn(specie.element, pw_read.structure.species, msg='Site {} not in read structure (sites: {})'.format(specie, pw_read.structure.species))
        # we construct this section dict manually because of the system parameters that are added from the struct
        sections_dict = {'control': {'calculation': 'scf', 'pseudo_dir': './'}, 'system': {'ecutwfc': 50.0, 'ibrav': 0, 'nat': 3, 'ntyp': 2}, 'electrons': {}, 'ions': {}, 'cell': {}}
        self.assertDictEqual(sections_dict, pw_read.sections, msg='Written and read parameter dictionaries are inconsistent')

class PWOuputTest(PymatgenTest):

    def setUp(self):
        self.pwout = PWOutput(os.path.join(test_dir, "Si.pwscf.out"))

    def test_properties(self):
        self.assertAlmostEqual(self.pwout.final_energy, -93.45259708)

    def test_get_celldm(self):
        self.assertAlmostEqual(self.pwout.get_celldm(1), 10.323)
        for i in range(2, 7):
            self.assertAlmostEqual(self.pwout.get_celldm(i), 0)

if __name__ == '__main__':
    unittest.main()
