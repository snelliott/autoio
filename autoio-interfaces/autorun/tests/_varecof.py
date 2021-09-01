""" tests varecof runners
"""

import os
import autorun
from ioformat import pathtools

PATH = os.path.dirname(os.path.realpath(__file__))

def test__frame_geoms():
    """ test autorun.frame__()
    """

    frag_geoms = autorun.varecof.frame_oriented_structure(
        script_str, run_dir, )

    for frag_geom, ref_geom in zip(frag_geoms, ref_geoms):
        assert automol.geom.almost_equal_dist_matrix(frag_geom, ref_geom)
#
#
# def test__compile_potentials():
#     """ test autorun.compile_potentials()
#     """
#
#     compile_potentials(run_dir, mep_distances, potentials,
#                        bnd_frm_idxs, fortran_compiler,
#                        dist_restrict_idxs=(),
#                        pot_labels=(),
#                        pot_file_names=(),
#                        spc_name='mol')
#
#     pot_str = pathtools.read_file(DAT_PATH, 'libcorrpot64.so')
#     ref_pot_str = pathtools.read_file(
#         DAT_PATH, 'libcorrpot64.so')
#     assert pot_str == ref_pot_str
#
#
#
if __name__ == '__main__':
#     test__run()
#     test__frame_geoms()
#     test__compile_potentials()
