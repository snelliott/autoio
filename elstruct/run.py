""" core run function
"""

from autorun import from_input_string
from autorun import dump_input
from autorun import dump_output
from autorun.ase import from_calc_dictionary as ase_from_dictionary

def direct(input_writer, script_str, run_dir, prog,
           geo, charge, mult, method, basis, **kwargs):
    """ Generates an input file for an electronic structure job and
        runs it directly.

        :param input_writer: elstruct writer module function for desired job
        :type input_writer: elstruct function
        :param script_str: string of bash script that contains
            execution instructions electronic structure job
        :type script_str: str
        :param run_dir: name of directory to run electronic structure job
        :type run_dir: str
        :param prog: electronic structure program to run
        :type prog: str
        :param geo: cartesian or z-matrix geometry
        :type geo: tuple
        :param charge: molecular charge
        :type charge: int
        :param mult: spin multiplicity
        :type mult: int
        :param method: electronic structure method
        :type method: str
        :returns: the input string, the output string, and the run directory
        :rtype: (str, str)
    """

    input_obj = input_writer(
        prog=prog,
        geo=geo, charge=charge, mult=mult, method=method, basis=basis,
        **kwargs)
    print(f"Input file for {prog}:\n{input_obj}\n")
    if isinstance(input_obj, str):
        output_strs = from_input_string(script_str, run_dir, input_obj)
        output_obj = output_strs[0]
        
    elif isinstance(input_obj, dict):
        output_obj = api(input_obj, script_str, run_dir)
        input_str = "" # figure out what to return later

    return input_obj, output_obj


def api(input_dct, script_str, run_dir):
    """ Generates an input file for an electronic structure job and
        runs it directly.

        :param input_writer: elstruct writer module function for desired job
        :type input_writer: elstruct function
        :param script_str: string of bash script that contains
            execution instructions electronic structure job
        :type script_str: str
        :param run_dir: name of directory to run electronic structure job
        :type run_dir: str
        :returns: the input string, the output string, and the run directory
        :rtype: (str, str)
    """
    dump_input(run_dir, input_dct)

    if 'ase' in script_str:
        output_dct = ase_from_dictionary(input_dct, script_str)

    dump_output(run_dir, output_dct)
    print(f"Output dictionary:\n{output_dct}\n")
    return output_dct
