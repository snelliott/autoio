""" MESS
"""

import mess_io.writer
from autorun._run import from_input_string


INPUT_NAME = 'mess.inp'
OUTPUT_NAMES = ('rate.out', 'mess.aux', 'mess.log')


# Specilialized runners
def torsions(script_str, run_dir, geo, hind_rot_str):
    """ Calculate the frequencies and ZPVES of the hindered rotors
        create a messpf input and run messpf to get tors_freqs and tors_zpes

        :param minthresh: threshold to assess if fits done to potential by MESS
            produce singularities
        :type minthresh: float
    """

    # Write the MESSPF input file
    input_str = mess_io.writer.messhr_inp_str(geo, hind_rot_str)

    # Run the direct function
    input_name = 'pf.inp'
    output_names = ('pf.out', 'pf.log', 'pf.aux')
    output_strs = direct(script_str, run_dir, input_str,
                         aux_dct=None,
                         input_name=input_name,
                         output_names=output_names)

    # Read the torsional freqs from output file
    out_str = output_strs[2]
    # tors_freqs = mess_io.reader.tors.first_point_harmonic_frequencies(out_str)

    # Read the torsional freqs and zpves from log file
    log_str = output_strs[1]
    tors_freqs = mess_io.reader.tors.analytic_frequencies(log_str)
    tors_zpes = mess_io.reader.tors.zero_point_vibrational_energies(
        log_str)

    # Read the fitted potential, report errors that are found
    # tors_pot_mins = mess_io.reader....()
    # for pot_min in tors_pot_mins:
    #     if pot_min < minthresh:
    #         print('WARNING: Potential fit from MESS has singularities')

    return tors_freqs, tors_zpes


def direct(script_str, run_dir, input_str, aux_dct=None,
           input_name=INPUT_NAME,
           output_names=OUTPUT_NAMES):
    """
        :param aux_dct: auxiliary input strings dict[name: string]
        :type aux_dct: dict[str: str]
        :param script_str: string of bash script that contains
            execution instructions electronic structure job
        :type script_str: str
        :param run_dir: name of directory to run electronic structure job
        :type run_dir: str
    """

    output_strs = from_input_string(
        script_str, run_dir, input_str,
        aux_dct=aux_dct,
        input_name=input_name,
        output_names=output_names)

    return output_strs
