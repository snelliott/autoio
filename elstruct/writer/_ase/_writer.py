""" ASE writer module """

import os
from ioformat import build_mako_str
import automol
import elstruct.option
import elstruct.par
from elstruct.writer import fill
from elstruct.writer._ase._par import OPTION_EVAL_DCT

PROG = elstruct.par.Program.ASE

def write_input(job_key, geo, charge, mult, method, basis, orb_restricted,
                # molecule options
                mol_options=(),
                # machine options
                memory=1, comment='', machine_options=(),
                # theory options
                scf_options=(), casscf_options=(), corr_options=(),
                # generic options
                gen_lines=None,
                # job options
                job_options=(), frozen_coordinates=(), saddle=False):
    """ Write an input file string for an electronic structure calculation
        by processing all of the information and using it to fill in
        a Mako template of the input file.

        :param job_key: job contained in the input file
        :type job_key: str
        :param geo: cartesian or z-matrix geometry
        :type geo: tuple
        :param charge: molecular charge
        :type charge: int
        :param mult: spin multiplicity
        :type mult: int
        :param method: electronic structure method
        :type method: str
        :param basis: basis set
        :type basis: str
        :param orb_restricted: parameter designating if restriced reference used
        :type orb_restricted: bool
        :param mol_options: options for the molecule block
        :type mol_options: tuple[str]
        ;param memory: memory in GB
        :type memory: int
        :param comment: a comment string to be placed at the top of the file
        :type comment: str
        :param machine_options: machine directives
            (num procs, num threads, etc.)
        :type machine_options: tuple[str]
        :param scf_options: scf method directives
        :type scf_options: tuple[str]
        :param casscf_options: casscf method directives
        :type casscf_options: tuple[str]
        :param corr_options: correlation method directives
        :type corr_options: tuple[str]
        :param job_options: geometry optimization routine directives
        :type job_options: tuple[str]
        :param frozen_coordinates: only with z-matrix geometries; list of
            coordinate names to freeze
        :type fozen_coordinates: tuple[str]
        :param saddle: parameter signifiying a saddle point calculation
        :type saddle: bool
        :param gen_lines: generic lines for the input file
        :type gen_lines: dict[idx:str]
    """
    calc_parameters = {}
    atom_parameters = to_ase_atom_object_parameters(geo)
    calc_parameters['atom_parameters'] = atom_parameters
    calc_parameters['job'] = job_key
    calc_parameters['charge'] = charge
    calc_parameters['multiplicity'] = mult
    calc_parameters['method'] = method
    calc_parameters['basis'] = basis
    calc_parameters['reference'] = 'rhf' if orb_restricted else 'uhf'

    print('calc_parameters:', calc_parameters)  # Debug print statement
    return calc_parameters 


def to_ase_atom_object_parameters(geo):
    atom_parameters = {}
    if automol.zmat.is_valid(geo):
        geo = automol.zmat.geometry(geo)
    atom_parameters['symbols'] = automol.geom.symbols(geo)
    atom_parameters['positions'] = automol.geom.coordinates(geo, angstrom=True)
    return atom_parameters



