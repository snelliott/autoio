""" test mess_io.reader.pes
"""

import os
import numpy
import mess_io
from ioformat import pathtools
from ioformat import remove_comment_lines
from autofile.io_ import read_file
import autoparse.pattern as app


PATH = os.path.dirname(os.path.realpath(__file__))
INP_PATH = os.path.join(PATH, 'data', 'inp')
INP_STR = pathtools.read_file(INP_PATH, 'mess.inp')


def test_pes():
    """ test mess_io.reader.pes
    """

    # Test reading with removing any fake wells
    energy_dct1, conn_lst1, _, _ = mess_io.reader.pes(
        input_string=INP_STR,
        read_fake=False)

    ref_energy_dct1 = {
        'P1': 0.0,
        'P2': 3.22,
        'B1': 13.23
    }
    ref_conn_lst1 = (
        ('P1', 'B1'),
        ('B1', 'P2')
    )

    ref_keys = tuple(ref_energy_dct1.keys())
    assert set(energy_dct1.keys()) == set(ref_keys)
    for key in ref_keys:
        assert numpy.isclose(energy_dct1[key], ref_energy_dct1[key])
    assert conn_lst1 == ref_conn_lst1

    # Test reading the entire PES with fake wells
    energy_dct2, conn_lst2, _, _ = mess_io.reader.pes(
        input_string=INP_STR,
        read_fake=True)

    ref_energy_dct2 = {
        'F1': -1.0,
        'F2': 2.22,
        'P1': 0.0,
        'P2': 3.22,
        'FRB1': 0.0,
        'FPB1': 3.22,
        'B1': 13.23
    }

    ref_conn_lst2 = (
        ('P1', 'FRB1'),
        ('FRB1', 'F1'),
        ('P2', 'FPB1'),
        ('FPB1', 'F2'),
        ('F1', 'B1'),
        ('B1', 'F2')
    )

    ref_keys = tuple(ref_energy_dct2.keys())
    assert set(energy_dct2.keys()) == set(ref_keys)
    for key in ref_keys:
        assert numpy.isclose(energy_dct2[key], ref_energy_dct2[key])
    assert conn_lst2 == ref_conn_lst2


def test_get_species():
    me_ped_inp = read_file(os.path.join(
        PATH, 'test_nB', 'C3H8_H', 'me_ktp_ped.inp'))
    me_ped_inp = remove_comment_lines(
        me_ped_inp, delim_pattern=app.escape('!'))
    me_ped_inp = remove_comment_lines(
        me_ped_inp, delim_pattern=app.escape('#'))
    species_blocks_ped = mess_io.reader.get_species(me_ped_inp)
    assert list(species_blocks_ped.keys()) == ['W0', 'RH', 'NC3H7', 'IC3H7']
    assert [len(i) for i in species_blocks_ped.values()] == [1, 2, 2, 2]


def test_find_barrier():
    conn_lst_dct = {'B0': ('W0', 'RH'), 'B1': ('W0', 'NC3H7'), 'B2': ('W0', 'IC3H7')}
    barrier_label1 = mess_io.reader.find_barrier(conn_lst_dct, 'RH', 'NC3H7')
    assert barrier_label1 == None
    barrier_label2 = mess_io.reader.find_barrier(conn_lst_dct, 'RH', 'W0')
    assert barrier_label2 == 'B0'
    barrier_label3 = mess_io.reader.find_barrier(conn_lst_dct, 'W0', 'NC3H7')
    assert barrier_label3 == 'B1'


def test_dct_species_fragments():
    # first extract species blocks
    me_ped_inp = read_file(os.path.join(
        PATH, 'test_nB', 'C3H8_H', 'me_ktp_ped.inp'))
    me_ped_inp = remove_comment_lines(
        me_ped_inp, delim_pattern=app.escape('!'))
    me_ped_inp = remove_comment_lines(
        me_ped_inp, delim_pattern=app.escape('#'))
    species_blocks_ped = mess_io.reader.get_species(me_ped_inp)
    dct_sp_fr = mess_io.reader.dct_species_fragments(species_blocks_ped)
    assert dct_sp_fr == {'W0': ['W0'], 'RH': ['C3H8', 'H'], 'NC3H7': [
        'CH3CH2CH2', 'H2'], 'IC3H7': ['CH3CHCH3', 'H2']}


if __name__ == '__main__':
    test_get_species()
    test_dct_species_fragments()
    test_find_barrier()

