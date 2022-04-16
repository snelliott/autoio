""" test mess_io.reader.hoten
"""

import os
import numpy as np
from ioformat import pathtools
import mess_io


PATH = os.path.dirname(os.path.realpath(__file__))
IPATH = os.path.join(PATH, 'data', 'inp')
OPATH = os.path.join(PATH, 'data', 'out')
HOT_INP_SGL = pathtools.read_file(IPATH, 'me_ktp_hoten_ch2o_oh.inp')
HOT_LOG_SGL = pathtools.read_file(OPATH, 'me_ktp_hoten_ch2o_oh.logf')
HOT_INP_DBL = pathtools.read_file(IPATH, 'me_ktp_hoten_c3h7.inp')
HOT_LOG_DBL = pathtools.read_file(OPATH, 'me_ktp_hoten_c3h7.logf')


def test_get_hot_species():
    """ test mess_io.read.get_hot_species
    """
    assert mess_io.reader.hoten.get_hot_species(HOT_INP_SGL) == {
        'W1': 0.0}
    assert mess_io.reader.hoten.get_hot_species(HOT_INP_DBL) == {
        'CH3CH2CH2': 3.19, 'CH3CHCH3': 0.0}


def test_extract_hot_branching():
    """ test mess_io.read.extract_hot_branching
    """
    hotspecies_en = {'W1': 0.0}
    species_lst = ('W1', 'CO+H')
    hoten_branch_dct = mess_io.reader.hoten.extract_hot_branching(
        HOT_LOG_SGL, hotspecies_en, species_lst)
    hoten = hoten_branch_dct['W1']

    assert np.allclose(hoten[3.16][1200].iloc[0].values,
                       np.array([0, 1]))
    assert np.allclose(hoten[3.16][1200].iloc[-1].values,
                       np.array([0.99764085, 0.00235915]))
    assert np.allclose(hoten[3.16][1200].loc[21.20].values,
                       np.array([5.939647e-05, 0.999941]))

    hotspecies_en = {'CH3CH2CH2': 3.19, 'CH3CHCH3': 0.0}
    species_lst = ('CH3CH2CH2', 'CH3CHCH3', 'C2H4+CH3', 'CH3CHCH2+H')
    hoten_branch_dct = mess_io.reader.hoten.extract_hot_branching(
        HOT_LOG_DBL, hotspecies_en, species_lst)
    hoten = hoten_branch_dct['CH3CHCH3']

    assert np.allclose(hoten[100][1800].iloc[0].values,
                       np.array([1.08041155e-06, 1.18044966e-04,
                                 4.65177195e-02, 9.53363155e-01]),
                       atol=1e-5)
    assert np.allclose(hoten[100][1800].iloc[-1].values,
                       np.array([4.28816038e-04, 9.99571184e-01,
                                 0.00000000e+00, 0.00000000e+00]),
                       atol=1e-5)
    assert np.allclose(hoten[100][1800].loc[53.3].values,
                       np.array([0.000473,  0.992253,  0.000000,    0.007275]),
                       atol=1e-5)

def test_extract_fne():

    dct_bf_tp_df = mess_io.reader.hoten.extract_fne(HOT_LOG_SGL)
    assert np.allclose((dct_bf_tp_df['W1'][0.1][1000]).values,
                       np.array([0.988, 0.012]))

    dct_bf_tp_df = mess_io.reader.hoten.extract_fne(HOT_LOG_DBL)
    assert np.allclose((dct_bf_tp_df['CH3CH2CH2'][0.1][1000]).values,
                       np.array([9.913503e-01, 6.542312e-07, 8.543019e-03, 1.060375e-04]))
    assert np.allclose((dct_bf_tp_df['CH3CHCH3'][0.1][1000]).values,
                       np.array([2.219028e-07, 9.995622e-01, 1.729243e-06, 4.358091e-04]))

if __name__ == '__main__':
    test_get_hot_species()
    test_extract_hot_branching()
    test_extract_fne()
