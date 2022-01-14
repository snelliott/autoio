""" helpers for importing and managing program modules
"""

import importlib
from elstruct import par
from elstruct import pclass


# Functions to import and call the appropriate writer function
def call_module_function(prog, function, *args, **kwargs):
    """ call the module implementation of a given function

        :param prog: the program
        :type prog: str
        :param function_template: a function with the desired signature
        :type function_template: function
    """

    def _rename_prog(prog):
        """ Rename a program if number does not match module name """
        if prog == 'molpro2021':
            prog = 'molpro2015'
        elif prog == 'gaussian16':
            prog = 'gaussian09'
        return prog

    assert prog in pclass.values(par.Program)
    assert prog in program_modules_with_function(function)

    name = f'_{_rename_prog(prog)}'
    module = importlib.import_module(f'elstruct.reader.{name:s}')
    reader = getattr(module, function)

    return reader(*args, **kwargs)


def program_modules_with_function(function):
    """
        :param function: a function with the desired signature
        :type function: function
    """

    progs = []
    for prog in pclass.values(par.Program):
        if function in READER_MODULE_DCT[prog]:
            progs.append(prog)

    return progs


# Information on what writers have been implemented
class Job():
    """ Names of electronic structure jobs to ne written
    """
    ENERGY = 'energy'
    GRADIENT = 'gradient'
    HESSIAN = 'hessian'
    HARM_FREQS = 'harmonic_frequencies'
    NORM_COORDS = 'normal_coordinates'
    IRC_PTS = 'irc_points'
    IRC_PATH = 'irc_path'
    OPT_GEO = 'opt_geometry'
    OPT_ZMA = 'opt_zmatrix'
    INP_ZMA = 'inp_zmatrix'
    VPT2 = 'vpt2'
    DIP_MOM = 'dipole_moment'
    POLAR = 'polarizability'
    EXIT_MSG = 'has_normal_exit_message'
    ERR_LST = 'error_list'
    SUCCESS_LST = 'success_list'
    ERR_MSG = 'has_error_message'
    CONV_MSG = 'check_convergence_messages'
    PROG_NAME = 'program_name'
    PROG_VERS = 'program_version'


# Dictionaries that dictate what writer/reader functionality
READER_MODULE_DCT = {
    par.Program.CFOUR2: (
        Job.ENERGY, Job.GRADIENT,
        Job.OPT_GEO, Job.OPT_ZMA,
        Job.EXIT_MSG, Job.ERR_LST, Job.SUCCESS_LST,
        Job.ERR_MSG, Job.CONV_MSG,
        Job.PROG_NAME, Job.PROG_VERS
    ),
    par.Program.GAUSSIAN09: (
        Job.ENERGY, Job.GRADIENT,
        Job.HESSIAN, Job.HARM_FREQS, Job.NORM_COORDS,
        Job.IRC_PTS, Job.IRC_PATH,
        Job.OPT_GEO, Job.OPT_ZMA, Job.INP_ZMA,
        Job.VPT2,
        Job.DIP_MOM, Job.POLAR,
        Job.EXIT_MSG, Job.ERR_LST, Job.SUCCESS_LST,
        Job.ERR_MSG, Job.CONV_MSG,
        Job.PROG_NAME, Job.PROG_VERS
    ),
    par.Program.GAUSSIAN16: (
        Job.ENERGY, Job.GRADIENT,
        Job.HESSIAN, Job.HARM_FREQS, Job.NORM_COORDS,
        Job.IRC_PTS, Job.IRC_PATH,
        Job.OPT_GEO, Job.OPT_ZMA,
        Job.VPT2,
        Job.DIP_MOM, Job.POLAR,
        Job.EXIT_MSG, Job.ERR_LST, Job.SUCCESS_LST,
        Job.ERR_MSG, Job.CONV_MSG,
        Job.PROG_NAME, Job.PROG_VERS
    ),
    par.Program.MOLPRO2015: (
        Job.ENERGY, Job.GRADIENT,
        Job.HESSIAN, Job.HARM_FREQS, Job.NORM_COORDS,
        Job.OPT_GEO, Job.OPT_ZMA,
        Job.EXIT_MSG, Job.ERR_LST,
        Job.ERR_MSG, Job.CONV_MSG,
        Job.PROG_NAME, Job.PROG_VERS
    ),
    par.Program.MOLPRO2021: (
        Job.ENERGY, Job.GRADIENT,
        Job.HESSIAN, Job.HARM_FREQS, Job.NORM_COORDS,
        Job.OPT_GEO, Job.OPT_ZMA,
        Job.EXIT_MSG, Job.ERR_LST,
        Job.ERR_MSG, Job.CONV_MSG,
        Job.PROG_NAME, Job.PROG_VERS
    ),
    par.Program.MRCC2018: (
        Job.ENERGY, Job.GRADIENT,
        Job.DIP_MOM,
        Job.EXIT_MSG, Job.ERR_LST,
        Job.ERR_MSG, Job.CONV_MSG,
        Job.PROG_NAME, Job.PROG_VERS
    ),
    par.Program.NWCHEM6: (
        Job.ENERGY, Job.GRADIENT,
        Job.OPT_GEO,
        Job.PROG_NAME, Job.PROG_VERS
    ),
    par.Program.ORCA4: (
        Job.ENERGY, Job.GRADIENT,
        Job.HESSIAN,
        Job.OPT_GEO,
        Job.DIP_MOM,
        Job.EXIT_MSG, Job.ERR_LST,
        Job.ERR_MSG, Job.CONV_MSG,
        Job.PROG_NAME, Job.PROG_VERS
    ),
    par.Program.PSI4: (
        Job.ENERGY, Job.GRADIENT,
        Job.HESSIAN, Job.HARM_FREQS,
        Job.IRC_PTS, Job.IRC_PATH,
        Job.OPT_GEO, Job.OPT_ZMA, Job.INP_ZMA,
        Job.DIP_MOM, Job.POLAR,
        Job.EXIT_MSG, Job.ERR_LST,
        Job.ERR_MSG, Job.CONV_MSG,
        Job.PROG_NAME, Job.PROG_VERS
    ),
    par.Program.QCHEM5: ()
}
