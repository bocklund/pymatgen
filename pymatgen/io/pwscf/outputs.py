from collections import defaultdict
from xml.etree import cElementTree as ET

from scipy.constants import physical_constants
from monty.re import regrep
from pymatgen import Lattice, Structure

bohr_to_m = physical_constants['Bohr radius'][0]
bohr_to_a = bohr_to_m*1e10

class PWOutput(object):

    patterns = {
        "energies": r'total energy\s+=\s+([\d\.\-]+)\sRy',
        "ecut": r'kinetic\-energy cutoff\s+=\s+([\d\.\-]+)\s+Ry',
        "lattice_type": r'bravais\-lattice index\s+=\s+(\d+)',
        "celldm1": r"celldm\(1\)=\s+([\d\.]+)\s",
        "celldm2": r"celldm\(2\)=\s+([\d\.]+)\s",
        "celldm3": r"celldm\(3\)=\s+([\d\.]+)\s",
        "celldm4": r"celldm\(4\)=\s+([\d\.]+)\s",
        "celldm5": r"celldm\(5\)=\s+([\d\.]+)\s",
        "celldm6": r"celldm\(6\)=\s+([\d\.]+)\s",
        "nkpts": r"number of k points=\s+([\d]+)"
    }

    def __init__(self, filename):
        self.filename = filename
        self.data = defaultdict(list)
        self.read_pattern(PWOutput.patterns)
        for k, v in self.data.items():
            if k == "energies":
                self.data[k] = [float(i[0][0]) for i in v]
            elif k in ["lattice_type", "nkpts"]:
                self.data[k] = int(v[0][0][0])
            else:
                self.data[k] = float(v[0][0][0])

    def read_pattern(self, patterns, reverse=False,
                     terminate_on_match=False, postprocess=str):
        """
        General pattern reading. Uses monty's regrep method. Takes the same
        arguments.

        Args:
            patterns (dict): A dict of patterns, e.g.,
                {"energy": r"energy\\(sigma->0\\)\\s+=\\s+([\\d\\-.]+)"}.
            reverse (bool): Read files in reverse. Defaults to false. Useful for
                large files, esp OUTCARs, especially when used with
                terminate_on_match.
            terminate_on_match (bool): Whether to terminate when there is at
                least one match in each key in pattern.
            postprocess (callable): A post processing function to convert all
                matches. Defaults to str, i.e., no change.

        Renders accessible:
            Any attribute in patterns. For example,
            {"energy": r"energy\\(sigma->0\\)\\s+=\s+([\\d\\-.]+)"} will set the
            value of self.data["energy"] = [[-1234], [-3453], ...], to the
            results from regex and postprocess. Note that the returned
            values are lists of lists, because you can grep multiple
            items on one line.
        """
        matches = regrep(self.filename, patterns, reverse=reverse,
                         terminate_on_match=terminate_on_match,
                         postprocess=postprocess)
        self.data.update(matches)

    def get_celldm(self, i):
        return self.data["celldm%d" % i]

    @property
    def final_energy(self):
        return self.data["energies"][-1]

    @property
    def lattice_type(self):
        return self.data["lattice_type"]


class PWData(PWOutput):
    """Open and provide an interface for parsing the XML
    """

    def __init__(self, xml_filename, output_filename):
        with open(xml_filename) as f:
            parsed_xml = ET.parse(f)
        self.root = parsed_xml.getroot()

        # properties that can be calculated
        self._structure = None
        # get a PWOutput from the output file
        super().__init__(output_filename)

    @staticmethod
    def _format_str(string):
        """Convert a string with excess spaces and line endings to a float
        """
        return [x for x in string.replace('\n', '').split(' ') if x != ''][0]

    @property
    def structure(self):
        if self._structure is not None:
            return self._structure
        else:
            # get lattice vectors as direct coordinates (units Bohr)
            cell = self.root.find('CELL')
            lattice_vectors = cell.find('DIRECT_LATTICE_VECTORS').getchildren()[
                              1:]
            lattice_vectors = [[float(val) * bohr_to_a for val in
                                v.text.replace('\n', '').split(' ') if
                                val != ''] for v in lattice_vectors]
            # get atomic positions as cartesian coordinates (units Bohr)
            atoms = []
            positions = []
            for el in self.root.find('IONS').getchildren():
                if 'ATOM.' in el.tag:
                    attribs = el.attrib
                    atoms.append(attribs['SPECIES'].replace(' ', ''))
                    position = [float(val) * bohr_to_a for val in
                                attribs['tau'].split(' ') if val != '']
                    positions.append(position)
            # create the structure:
            lattice = Lattice(lattice_vectors)
            self._structure = Structure(lattice, atoms, positions,
                                        coords_are_cartesian=True)
            return self._structure

    def as_dict(self):
        natoms = len(self.structure.species)
        return {
            'formula_pretty': self.structure.composition.reduced_formula,
            'formula': self.structure.formula,
            'natoms': natoms,
            'kpoints': {
                'grid': tuple(int(v) for _, v in
                              self.root.find('BRILLOUIN_ZONE').find(
                                  'MONKHORST_PACK_GRID').attrib.items()),
                'shift': tuple(int(v) for _, v in
                               self.root.find('BRILLOUIN_ZONE').find(
                                   'MONKHORST_PACK_OFFSET').attrib.items()),
                'nkpoints': int(self._format_str(
                    self.root.find('BRILLOUIN_ZONE').find(
                        'NUMBER_OF_K-POINTS').text))
            },
            'output': {
                'wfc_cutoff': float(self._format_str(
                    self.root.find('PLANE_WAVES').find('WFC_CUTOFF').text)),
                'rho_cutoff': float(self._format_str(
                    self.root.find('PLANE_WAVES').find('RHO_CUTOFF').text)),
                'number_of_bands': int(self._format_str(
                    self.root.find('BAND_STRUCTURE_INFO').find(
                        'NUMBER_OF_BANDS').text)),
                'number_of_electrons': float(self._format_str(
                    self.root.find('BAND_STRUCTURE_INFO').find(
                        'NUMBER_OF_ELECTRONS').text)),
                'fermi_energy': float(self._format_str(
                    self.root.find('BAND_STRUCTURE_INFO').find(
                        'FERMI_ENERGY').text)),
                'units_energy': self._format_str(
                    self.root.find('BAND_STRUCTURE_INFO').find(
                        'UNITS_FOR_ENERGIES').attrib['UNITS']),
                'structure': self.structure.as_dict(),
                'density': self.structure.density,
                'energy': self.final_energy / 2.0,  # to Hartrees
                'energy_per_atom': self.final_energy / 2.0 / natoms,
            }
        }

    @property
    def energy_per_atom(self):
        return self.final_energy / 2.0 / len(self.structure.species)
