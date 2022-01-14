""" test the writing of the energy transfer section
"""

import os
from ioformat import pathtools
import mess_io.writer


PATH = os.path.dirname(os.path.realpath(__file__))
INP_PATH = os.path.join(PATH, 'data', 'inp')


def test__energy_trans_writer():
    """ tests writing the section to a file
    """

    edown_str = mess_io.writer.energy_down(
        exp_factor=150.0,
        exp_power=50.0,
        exp_cutoff=80.0
    )
    assert edown_str == pathtools.read_file(
        INP_PATH, 'etrans_edown.inp')

    collid_str = mess_io.writer.collision_frequency(
        eps1=100.0,
        eps2=200.0,
        sig1=10.0,
        sig2=20.0,
        mass1=15.0,
        mass2=25.0)
    assert collid_str == pathtools.read_file(
        INP_PATH, 'etrans_collid.inp')

    glob_etrans_str = mess_io.writer.global_energy_transfer_input(
        edown_str, collid_str)

    assert glob_etrans_str == pathtools.read_file(
        INP_PATH, 'glob_etrans.inp')
