import os

from monty.json import MSONable

from pymatgen.io.pwscf.inputs import PWInput

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



# TODO: calculate optimal Kpoint density
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
        filename = os.path.join(path, self.structure.composition.reduced_formula + '.in')
        self.pwinput.write_file(filename)

    def as_dict(self, verbosity=2):
        d = MSONable.as_dict(self)
        if verbosity == 1:
            d.pop("structure", None)
        return d


class PWStaticSet(PWInputSet):
    """Return a PWInput object for a typical PWscf static calculation
    """
    def __init__(self, structure, pseudo=None, pseudo_dir=None, user_input_settings=None, user_kpoints_settings=None):
        pw_static_dict = {
            'control': {
                'calculation': 'scf',
                'prefix': structure.composition.reduced_formula,
                'outdir': './',
                'etot_conv_thr': 10 ** -6,
                'forc_conv_thr': 10 ** -5,
                'pseudo_dir': pseudo_dir
            },
            'system': {
                'ecutwfc': 130,
                # TODO: set this automatically. Determine the cutoffs independently
                'ecutrho': 6 * 130,
                'ntyp': len(structure.types_of_specie),
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
                'calculation': 'relax',
                'prefix': structure.composition.reduced_formula,
                'outdir': './',
                'nstep': 100,
                'etot_conv_thr': 10 ** -6,
                'forc_conv_thr': 10 ** -5,
                'pseudo_dir': pseudo_dir
            },
            'system': {
                'ecutwfc': 200,
            # TODO: set this automatically. Determine the cutoffs independently
                'ecutrho': 6 * 200,
                'ntyp': len(structure.types_of_specie),
            },
            'kpoints_shift': (0, 0, 0)
        }
        super(PWRelaxSet, self).__init__(structure, pseudo=pseudo, pw_input_dict=pw_relax_dict, pseudo_dir=pseudo_dir, user_input_settings=user_input_settings, user_kpoints_settings=user_kpoints_settings)
