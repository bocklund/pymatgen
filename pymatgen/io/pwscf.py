# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

import os
import re
import six
import xml.etree.cElementTree as ET

from monty.io import zopen
from monty.json import MSONable
from monty.re import regrep
from collections import defaultdict
from scipy.constants import physical_constants

from pymatgen.core.periodic_table import Element
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.util.io_utils import clean_lines

"""
This module implements input and output processing from PWSCF.
"""

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Virtual Lab"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "ongsp@ucsd.edu"
__date__ = "3/27/15"

bohr_to_m = physical_constants['Bohr radius'][0]
bohr_to_a = bohr_to_m*1e10


SSSP_accurate = {
'Al': 'Al.pbe-n-kjpaw_psl.1.0.0.UPF',
'Os': 'Os.pbe-spfn-rrkjus_psl.1.0.0.UPF',
'Ar': 'Ar.pbe-n-rrkjus_psl.1.0.0.UPF',
'P': 'P.pbe-n-rrkjus_psl.1.0.0.UPF',
'As': 'As.pbe-n-rrkjus_psl.0.2.UPF',
'Pb': 'Pb.pbe-dn-kjpaw_psl.0.2.2.UPF',
'Au': 'Au_ONCV_PBE-1.0.upf',
'Pd': 'Pd.pbe-spn-kjpaw_psl.1.0.0.UPF',
'B': 'B.pbe-n-kjpaw_psl.0.1.UPF',
'Pm': 'Pm.GGA-PBE-paw-v1.0.UPF',
'Ba': 'Ba_ONCV_PBE-1.0.upf',
'Po': 'Po.pbe-dn-rrkjus_psl.1.0.0.UPF',
'Be': 'Be_ONCV_PBE-1.0.upf',
'Pr': 'Pr.GGA-PBE-paw-v1.0.UPF',
'Bi': 'Bi.pbe-dn-kjpaw_psl.0.2.2.UPF',
'Pt': 'Pt.pbe-spfn-rrkjus_psl.1.0.0.UPF',
'C': 'C_pbe_v1.2.uspp.F.UPF',
'Ca': 'Ca_pbe_v1.uspp.F.UPF',
'Rb': 'Rb_ONCV_PBE-1.0.upf',
'Cd': 'Cd.pbe-dn-rrkjus_psl.0.3.1.UPF',
'Re': 'Re_pbe_v1.2.uspp.F.UPF',
'Ce': 'Ce.GGA-PBE-paw-v1.0.UPF',
'Rh': 'Rh.pbe-spn-kjpaw_psl.1.0.0.UPF',
'Cl': 'Cl.pbe-n-rrkjus_psl.1.0.0.UPF',
'Rn': 'Rn.pbe-dn-rrkjus_psl.1.0.0.UPF',
'Co': 'Co_pbe_v1.2.uspp.F.UPF',
'Ru': 'Ru_ONCV_PBE-1.0.upf',
'Cs': 'Cs_pbe_v1.uspp.F.UPF',
'S': 'S_pbe_v1.2.uspp.F.UPF',
'Cu': 'Cu_pbe_v1.2.uspp.F.UPF',
'Sc': 'Sc_pbe_v1.uspp.F.UPF',
'Dy': 'Dy.GGA-PBE-paw-v1.0.UPF',
'Se': 'Se_pbe_v1.uspp.F.UPF',
'Er': 'Er.GGA-PBE-paw-v1.0.UPF',
'Si': 'Si.pbe-n-rrkjus_psl.1.0.0.UPF',
'Eu': 'Eu.GGA-PBE-paw-v1.0.UPF',
'Sm': 'Sm.GGA-PBE-paw-v1.0.UPF',
'Fe': 'Fe.pbe-spn-kjpaw_psl.0.2.1.UPF',
'Sn': 'Sn_pbe_v1.uspp.F.UPF',
'Ga': 'Ga.pbe-dn-kjpaw_psl.1.0.0.UPF',
'Sr': 'Sr.pbe-spn-rrkjus_psl.1.0.0.UPF',
'Gd': 'Gd.GGA-PBE-paw-v1.0.UPF',
'Ta': 'Ta.pbe-spfn-rrkjus_psl.1.0.0.UPF',
'Ge': 'Ge.pbe-dn-kjpaw_psl.1.0.0.UPF',
'Tb': 'Tb.GGA-PBE-paw-v1.0.UPF',
'H': 'H.pbe-rrkjus_psl.0.1.UPF',
'Tc': 'Tc_ONCV_PBE-1.0.upf',
'He': 'He_ONCV_PBE-1.0.upf',
'Te': 'Te_pbe_v1.uspp.F.UPF',
'Hf': 'Hf.pbe-spdfn-kjpaw_psl.1.0.0.UPF',
'Tl': 'Tl.pbe-dn-rrkjus_psl.1.0.0.UPF',
'Hg': 'Hg_pbe_v1.uspp.F.UPF',
'Tm': 'Tm.GGA-PBE-paw-v1.0.UPF',
'Ho': 'Ho.GGA-PBE-paw-v1.0.UPF',
'V': 'V_pbe_v1.uspp.F.UPF',
'I': 'I_pbe_v1.uspp.F.UPF',
'W': 'W_pbe_v1.2.uspp.F.UPF',
'In': 'In.pbe-dn-rrkjus_psl.0.2.2.UPF',
'Xe': 'Xe.pbe-dn-rrkjus_psl.1.0.0.UPF',
'Ir': 'Ir_pbe_v1.2.uspp.F.UPF',
'Y': 'Y_pbe_v1.uspp.F.UPF',
'K': 'K.pbe-spn-rrkjus_psl.1.0.0.UPF',
'Yb': 'Yb.GGA-PBE-paw-v1.0.UPF',
'Kr': 'Kr.pbe-n-rrkjus_psl.0.2.3.UPF',
'Zn': 'Zn_pbe_v1.uspp.F.UPF',
'La': 'La.GGA-PBE-paw-v1.0.UPF',
'Zr': 'Zr_pbe_v1.uspp.F.UPF',
'Lu': 'Lu.GGA-PBE-paw-v1.0.UPF',
'Ag': 'ag_pbe_v1.4.uspp.F.UPF',
'Mn': 'Mn.pbe-spn-kjpaw_psl.0.3.1.UPF',
'Br': 'br_pbe_v1.4.uspp.F.UPF',
'Mo': 'Mo_ONCV_PBE-1.0.upf',
'Cr': 'cr_pbe_v1.5.uspp.F.UPF',
'N': 'N.pbe.theos.UPF',
'F': 'f_pbe_v1.4.uspp.F.UPF',
'Na': 'Na_pbe_v1.uspp.F.UPF',
'Li': 'li_pbe_v1.4.uspp.F.UPF',
'Nb': 'Nb.pbe-spn-kjpaw_psl.0.3.0.UPF',
'Mg': 'mg_pbe_v1.4.uspp.F.UPF',
'Nd': 'Nd.GGA-PBE-paw-v1.0.UPF',
'Ni': 'ni_pbe_v1.4.uspp.F.UPF',
'Ne': 'Ne.pbe-n-kjpaw_psl.1.0.0.UPF',
'Sb': 'sb_pbe_v1.4.uspp.F.UPF',
'O': 'O.pbe-n-kjpaw_psl.0.1.UPF',
'Ti': 'ti_pbe_v1.4.uspp.F.UPF',
}

class PWInput(object):
    """
    Base input file class. Right now, only supports no symmetry and is
    very basic.
    """

    def __init__(self, structure, pseudo=None, control=None, system=None,
                 electrons=None, ions=None, cell=None, kpoints_mode="automatic",
                 kpoints_grid=(1, 1, 1),kpoints_shift=(0, 0, 0)):
        """
        Initializes a PWSCF input file.

        Args:
            structure (Structure): Input structure. For spin-polarized calculation,
                properties (e.g. {"starting_magnetization": -0.5, 
                "pseudo": "Mn.pbe-sp-van.UPF"}) on each site is needed instead of 
                pseudo (dict).
            pseudo (dict): A dict of the pseudopotentials to use. Default to None.
            control (dict): Control parameters. Refer to official PWSCF doc
                on supported parameters. Default to {"calculation": "scf"}
            system (dict): System parameters. Refer to official PWSCF doc
                on supported parameters. Default to None, which means {}.
            electrons (dict): Electron parameters. Refer to official PWSCF doc
                on supported parameters. Default to None, which means {}.
            ions (dict): Ions parameters. Refer to official PWSCF doc
                on supported parameters. Default to None, which means {}.
            cell (dict): Cell parameters. Refer to official PWSCF doc
                on supported parameters. Default to None, which means {}.
            kpoints_mode (str): Kpoints generation mode. Default to automatic.
            kpoints_grid (sequence): The kpoint grid. Default to (1, 1, 1).
            kpoints_shift (sequence): The shift for the kpoints. Defaults to
                (0, 0, 0).
        """
        self.structure = structure
        sections = {}
        sections["control"] = control or {"calculation": "scf"}
        sections["system"] = system or {}
        sections["electrons"] = electrons or {}
        sections["ions"] = ions or {}
        sections["cell"] = cell or {}
        if pseudo == None:
            for site in structure:
                try:
                    site.properties['pseudo']
                except KeyError:
                    raise PWInputError("Missing %s in pseudo specification!" 
                                       % site)
        else:
            for species in self.structure.composition.keys():
                if species.symbol not in pseudo:
                    raise PWInputError("Missing %s in pseudo specification!" 
                                       % species.symbol)
        self.pseudo = pseudo

        self.sections = sections
        self.kpoints_mode = kpoints_mode
        self.kpoints_grid = kpoints_grid
        self.kpoints_shift = kpoints_shift

    def __str__(self):
        out = []
        site_descriptions = {}

        if self.pseudo != None:
            site_descriptions = self.pseudo
        else:
            c = 1
            for site in self.structure:
                name = None
                for k, v in site_descriptions.items():
                    if site.properties == v:
                        name = k

                if name == None:
                    name = site.specie.symbol+str(c)
                    site_descriptions[name] = site.properties
                    c += 1

        def to_str(v):
            if isinstance(v, six.string_types):
                return "'%s'" % v
            elif isinstance(v, float):
                return "%s" % str(v).replace("e", "d")
            elif isinstance(v, bool):
                if v:
                    return ".TRUE."
                else:
                    return ".FALSE."
            return v

        for k1 in ["control", "system", "electrons", "ions", "cell"]:
            v1 = self.sections[k1]
            out.append("&%s" % k1.upper())
            sub = []
            for k2 in sorted(v1.keys()):
                if isinstance(v1[k2], list):
                    n = 1
                    for l in v1[k2][:len(site_descriptions)]:
                        sub.append("  %s(%d) = %s" % (k2, n, to_str(v1[k2][n-1])))
                        n += 1
                else:
                    sub.append("  %s = %s" % (k2, to_str(v1[k2])))
            if k1 == "system":
                if 'ibrav' not in self.sections[k1]:
                    sub.append("  ibrav = 0")
                if 'nat' not in self.sections[k1]:
                    sub.append("  nat = %d" % len(self.structure))
                if 'ntyp' not in self.sections[k1]:
                    sub.append("  ntyp = %d" % len(site_descriptions))
            sub.append("/")
            out.append(",\n".join(sub))

        out.append("ATOMIC_SPECIES")
        for k, v in sorted(site_descriptions.items(), key=lambda i: i[0]):
            e = re.match(r"[A-Z][a-z]?", k).group(0)
            if self.pseudo is not None:
                p = v
            else:
                p = v['pseudo']
            out.append("  %s  %.4f %s" % (k, Element(e).atomic_mass, p)) 

        out.append("ATOMIC_POSITIONS crystal")
        if self.pseudo is not None:
            for site in self.structure:
                out.append("  %s %.6f %.6f %.6f" % (site.specie.symbol, site.a,
                                                    site.b, site.c))
        else:
            for site in self.structure:
                name = None
                for k, v in sorted(site_descriptions.items(),
                                   key=lambda i: i[0]):
                    if v == site.properties:
                        name = k
                out.append("  %s %.6f %.6f %.6f" % (name, site.a, site.b, site.c))

        out.append("K_POINTS %s" % self.kpoints_mode)
        kpt_str = ["%s" % i for i in self.kpoints_grid]
        kpt_str.extend(["%s" % i for i in self.kpoints_shift])
        out.append("  %s" % " ".join(kpt_str))
        out.append("CELL_PARAMETERS angstrom")
        for vec in self.structure.lattice.matrix:
            out.append("  %f %f %f" % (vec[0], vec[1], vec[2]))
        return "\n".join(out)

    def write_file(self, filename):
        """
        Write the PWSCF input file.

        Args:
            filename (str): The string filename to output to.
        """
        with open(filename, "w") as f:
            f.write(self.__str__())

    @staticmethod
    def from_file(filename):
        """
        Reads an PWInput object from a file.

        Args:
            filename (str): Filename for file

        Returns:
            PWInput object
        """
        with zopen(filename, "rt") as f:
            return PWInput.from_string(f.read())

    @staticmethod
    def from_string(string):
        """
        Reads an PWInput object from a string.

        Args:
            string (str): PWInput string

        Returns:
            PWInput object
        """
        lines = list(clean_lines(string.splitlines()))

        def input_mode(line):
            if line[0] == "&":
                return ("sections", line[1:].lower())
            elif "ATOMIC_SPECIES" in line:
                return ("pseudo", )
            elif "K_POINTS" in line:
                return ("kpoints", line.split("{")[1][:-1])
            elif "CELL_PARAMETERS" in line or "ATOMIC_POSITIONS" in line:
                return ("structure", line.split("{")[1][:-1])
            elif line == "/":
                return None
            else:
                return mode

        sections = {"control": {}, "system": {}, "electrons": {}, 
                    "ions": {}, "cell":{}}
        pseudo = {}
        pseudo_index = 0
        lattice = []
        species = []
        coords = []
        structure = None
        site_properties = {"pseudo":[]}
        mode = None
        for line in lines:
            mode = input_mode(line)
            if mode == None:
                pass
            elif mode[0] == "sections":
                section = mode[1]
                m = re.match(r'(\w+)\(?(\d*?)\)?\s*=\s*(.*)', line)
                if m:
                    key = m.group(1).strip()
                    key_ = m.group(2).strip()
                    val = m.group(3).strip()
                    if key_ != "":
                        if sections[section].get(key, None) == None:
                            val_ = [0.0]*20 # MAX NTYP DEFINITION
                            val_[int(key_)-1] = PWInput.proc_val(key, val)
                            sections[section][key] = val_

                            site_properties[key] = []
                        else:
                            sections[section][key][int(key_)-1] = PWInput.proc_val(key, val) 
                    else:
                        sections[section][key] = PWInput.proc_val(key, val)

            elif mode[0] == "pseudo":
                m = re.match(r'(\w+)\s+(\d*.\d*)\s+(.*)', line)
                if m:
                    pseudo[m.group(1).strip()] = {}
                    pseudo[m.group(1).strip()]["index"] = pseudo_index
                    pseudo[m.group(1).strip()]["pseudopot"] = m.group(3).strip()
                    pseudo_index += 1
            elif mode[0] == "kpoints":
                m = re.match(r'(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)', line)
                if m:
                    kpoints_grid = (int(m.group(1)), int(m.group(2)), int(m.group(3)))
                    kpoints_shift = (int(m.group(4)), int(m.group(5)), int(m.group(6)))
                else:
                    kpoints_mode = mode[1]
            elif mode[0] == "structure":
                m_l = re.match(r'(-?\d+\.?\d*)\s+(-?\d+\.?\d*)\s+(-?\d+\.?\d*)', line)
                m_p = re.match(r'(\w+)\s+(-?\d+\.\d*)\s+(-?\d+\.?\d*)\s+(-?\d+\.?\d*)', line)
                if m_l:
                    lattice += [ float(m_l.group(1)), float(m_l.group(2)), float(m_l.group(3)) ]
                elif m_p:
                    site_properties["pseudo"].append(pseudo[m_p.group(1)]["pseudopot"])
                    species += [pseudo[m_p.group(1)]["pseudopot"].split(".")[0]]
                    coords += [[float(m_p.group(2)), float(m_p.group(3)), float(m_p.group(4))]]

                    for k, v in site_properties.items():
                        if k != "pseudo":
                            site_properties[k].append(sections['system'][k][pseudo[m_p.group(1)]["index"]])
                if mode[1] == "angstrom":
                    coords_are_cartesian = True
                elif mode[1] == "crystal":
                    coords_are_cartesian = False

        structure = Structure(Lattice(lattice), species, coords, 
                              coords_are_cartesian=coords_are_cartesian,
                              site_properties=site_properties)
        return PWInput(structure=structure, control=sections["control"],
                       system=sections["system"], electrons=sections["electrons"], 
                       ions=sections["ions"], cell=sections["cell"], kpoints_mode=kpoints_mode,
                       kpoints_grid=kpoints_grid, kpoints_shift=kpoints_shift)

    def proc_val(key, val):
        """
        Static helper method to convert PWINPUT parameters to proper type, e.g.,
        integers, floats, etc.

        Args:
            key: PWINPUT parameter key
            val: Actual value of PWINPUT parameter.
        """
        float_keys = ('etot_conv_thr','forc_conv_thr','conv_thr','Hubbard_U','Hubbard_J0','defauss',
                      'starting_magnetization',)

        int_keys = ('nstep','iprint','nberrycyc','gdir','nppstr','ibrav','nat','ntyp','nbnd','nr1',
                    'nr2','nr3','nr1s','nr2s','nr3s','nspin','nqx1','nqx2','nqx3','lda_plus_u_kind',
                    'edir','report','esm_nfit','space_group','origin_choice','electron_maxstep',
                    'mixing_ndim','mixing_fixed_ns','ortho_para','diago_cg_maxiter','diago_david_ndim',
                    'nraise','bfgs_ndim','if_pos','nks','nk1','nk2','nk3','sk1','sk2','sk3','nconstr')

        bool_keys = ('wf_collect','tstress','tprnfor','lkpoint_dir','tefield','dipfield','lelfield',
                     'lorbm','lberry','lfcpopt','monopole','nosym','nosym_evc','noinv','no_t_rev',
                     'force_symmorphic','use_all_frac','one_atom_occupations','starting_spin_angle',
                     'noncolin','x_gamma_extrapolation','lda_plus_u','lspinorb','london',
                     'ts_vdw_isolated','xdm','uniqueb','rhombohedral','realxz','block',
                     'scf_must_converge','adaptive_thr','diago_full_acc','tqr','remove_rigid_rot',
                     'refold_pos')

        def smart_int_or_float(numstr):
            if numstr.find(".") != -1 or numstr.lower().find("e") != -1:
                return float(numstr)
            else:
                return int(numstr)

        try:
            if key in bool_keys:
                if val.lower() == ".true.":
                    return True
                elif val.lower() == ".false.":
                    return False
                else:
                    raise ValueError(key + " should be a boolean type!")

            if key in float_keys:
                return float(re.search(r"^-?\d*\.?\d*d?-?\d*", val.lower()).group(0).replace("d", "e"))

            if key in int_keys:
                return int(re.match(r"^-?[0-9]+", val).group(0))

        except ValueError:
            pass

        try:
            val = val.replace("d","e")
            return smart_int_or_float(val)
        except ValueError:
            pass

        if "true" in val.lower():
            return True
        if "false" in val.lower():
            return False

        m = re.match(r"^[\"|'](.+)[\"|']$", val)
        if m:
            return m.group(1)



class PWInputError(BaseException):
    pass


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

class PWInputSet(MSONable):
    """
    Base class for input sets in PWscf.

    Except for the kpoints, which are currently not implemented (planned to be
    implemented as in VaspInputSet objects), all subclasses need to do are to
    pass the relevant pw_input_dict and reasonable defaults are supplied by this
    class.
    """
    def __init__(self, structure, pseudo=None, pw_input_dict=None, pseudo_dir=None, user_input_settings=None, user_kpoints_settings=None):
        """
        Create a PWInput object from the passed structure and pw_input_dict
        Args:
            structure (Structure): pymatgen Structure object
            pseudo (dict): Dictionary of pseudopotential files located in
                pseudo_dir formatted as {'Si': 'si.pbe-paw.upf'} for the element
                and the name of the pseuopotential.
            pw_input_dict (dict): Dictionary of input settings corresponding to
                the namelists defined as in the PWInput object
            pseudo_dir (str): String of the full path to the directory containing pseudopotentials.
            user_kpoints_settings (dict): Not yet implemented. Will be like VaspInputSets.
        """
        pw_input_dict = pw_input_dict or {}
        pseudo_dir = pseudo_dir or os.environ.get('PSEUDO_DIR') or os.path.expanduser('~/pseudo')
        # pass only the correct psedopotentials
        pseudo = pseudo or SSSP_accurate
        if user_input_settings:
            for key, val in user_input_settings.items():
                pw_input_dict[key].update(val)
        if user_kpoints_settings:
            raise NotImplementedError()
        else:
            # automatic kpoints mode and grid
            # TODO: set to 'gamma' if hexagonal structure
            pw_input_dict['kpoints_mode'] = 'automatic'
            pw_input_dict['kpoints_grid'] = (11, 11, 11)
        self.user_kpoints_settings = user_kpoints_settings
        self.pseudo = {str(sp): pseudo[str(sp)] for sp in structure.species}
        self.structure = structure
        self.pseudo_dir = pseudo_dir
        self.pw_input_dict = pw_input_dict
        self.user_input_settings = user_input_settings
        self.pwinput = PWInput(structure, pseudo=self.pseudo, **self.pw_input_dict)

    def write_input(self, path):
        filename = self.structure.composition.reduced_formula + '.in'
        self.pwinput.write_file(filename)

    def as_dict(self, verbosity=2):
        d = MSONable.as_dict(self)
        if verbosity == 1:
            d.pop("structure", None)
        return d

# TODO: calculate optimal Kpoint density
class PWStaticSet(PWInputSet):
    """Return a PWInput object for a typical PWscf static calculation
    """
    def __init__(self, structure, pseudo=None, pseudo_dir=None, user_input_settings=None, user_kpoints_settings=None):
        pw_static_dict = {
            'control': {
                'calculation': 'scf',  # do a static calculation
                'prefix': structure.composition.reduced_formula,
                'outdir': './',
                'etot_conv_thr': 10 ** -6,
            # energy convergence criteria for relaxation. Default 10^-4
                'forc_conv_thr': 10 ** -5,
            # force convergence criteria for relaxation Default 10^-3
                'pseudo_dir': pseudo_dir
            },
            'system': {
                # the highest non-noble gas recommended cutoff is Hf: 120 Ry, follwed by Fe: 95 Ry. We'll do 130 here and relax is 1.5* for safety
                'ecutwfc': 130,
                # TODO: set this automatically. Determine the cutoffs independently
                'ecutrho': 6 * 130,
                'ntyp': len(structure.types_of_specie),
                # default is 4*ecutwfc. Manual says should be higher for PAW, esp. if vacuum or non-linear core
            },
            'kpoints_shift': (0, 0, 0)
        }
        super(PWStaticSet, self).__init__(structure, pseudo=pseudo, pw_input_dict=pw_static_dict, pseudo_dir=pseudo_dir, user_input_settings=user_input_settings, user_kpoints_settings=user_kpoints_settings)

class PWRelaxSet(PWInputSet):
    """Return a PWInput object for a typical PWscf relaxation
    """
    def __init__(self, structure, pseudo=None, pseudo_dir=None, user_input_settings=None, user_kpoints_settings=None):
        pw_relax_dict = {
            'control': {
                'calculation': 'relax',  # do a variable cell relaxation
                'prefix': structure.composition.reduced_formula,
                'outdir': './',
                'nstep': 100,  # take 100 ionic steps
                'etot_conv_thr': 10 ** -6,
            # energy convergence criteria for relaxation. Default 10^-4
                'forc_conv_thr': 10 ** -5,
            # force convergence criteria for relaxation Default 10^-3
                'pseudo_dir': pseudo_dir
            },
            'system': {
                'ecutwfc': 200,
            # TODO: set this automatically. Determine the cutoffs independently
                'ecutrho': 6 * 200,
                'ntyp': len(structure.types_of_specie),
            # default is 4*ecutwfc. Manual says should be higher for PAW, esp. if vacuum or non-linear core
            },
            'kpoints_shift': (0, 0, 0)
        }
        super(PWRelaxSet, self).__init__(structure, pseudo=pseudo, pw_input_dict=pw_relax_dict, pseudo_dir=pseudo_dir, user_input_settings=user_input_settings, user_kpoints_settings=user_kpoints_settings)


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


