""" Calculate Total Energy Distribution and related quantities
"""

import intder_io
from autorun._run import from_input_string


# Default names of input and output files
INPUT_NAME = 'intder.inp'
OUTPUT_NAMES = ('intder.out',)


# Specialized Runners
def ted_zmatrix_coordinates(script_str, run_dir,
                            geo, zma, hess, mode_idx,
                            input_name=INPUT_NAME,
                            output_names=OUTPUT_NAMES):
    """ Calculate internal coordinate rep of mode from TED

        Note that geo and zma atom ordering must match.
        Currently, calculation does not support dummy atoms/linear segments
    """

    # Run INTDER
    output_strs = direct(script_str, run_dir, zma, geo, hess,
                         input_name=input_name, output_names=output_names)
    output_str = output_strs[0]

    # Get modes from TED output
    if output_strs is not None:
        # Read the internal coordinates and TED assignments from output
        intl_coords = intder_io.reader.internal_coordinates(output_str)
        ted_assign = intder_io.reader.ted_assignments(output_str)

        # Get the Z-Matrix names
        ted_zmat_names = intder_io.ted_zmatrix_coordinates(
            zma, mode_idx, intl_coords, ted_assign)
    else:
        ted_zmat_names = None

    return ted_zmat_names


# Generalized Runner
def direct(script_str, run_dir, zma, geo, hess,
           input_name=INPUT_NAME, output_names=OUTPUT_NAMES):
    """ Generates an input file for a ProjRot job, runs it directly, and
        obtains all of the possible output file strings
    """

    input_str = intder_io.writer.input_file(geo, zma)
    aux_dct = {'file15': intder_io.writer.cartesian_hessian_file(hess)}

    output_strs = from_input_string(
        script_str, run_dir, input_str,
        aux_dct=aux_dct,
        input_name=input_name,
        output_names=output_names)

    return output_strs
