""" test mess_io.reader.hoten
"""

import os
import numpy as np
import mess_io
from autofile.io_ import read_file

PATH = os.path.dirname(os.path.realpath(__file__))
INP_PATH_SINGLE = os.path.join(PATH, 'test_nB', 'CH2O_OH_ORIG')
HOT_INP_SINGLE = read_file(os.path.join(INP_PATH_SINGLE, 'me_ktp_hoten.inp'))
HOT_OUT_SINGLE = read_file(os.path.join(INP_PATH_SINGLE, 'me_ktp_hoten.log'))
INP_PATH_DOUBLE = os.path.join(PATH, 'test_nB', 'C3H7')
HOT_INP_DOUBLE = read_file(os.path.join(INP_PATH_DOUBLE, 'me_ktp_hoten.inp'))
HOT_OUT_DOUBLE = read_file(os.path.join(INP_PATH_DOUBLE, 'me_ktp_hoten.log'))


def test_get_hot_names():
    """ test mess_io.read.get_hot_names
    """
    assert mess_io.reader.hoten.get_hot_names(HOT_INP_SINGLE) == ['W1']
    assert mess_io.reader.hoten.get_hot_names(HOT_INP_DOUBLE) == ['CH3CH2CH2', 'CH3CHCH3']

def test_extract_hot_branching():
    """ test mess_io.read.extract_hot_branching
    """
    hotspecies_lst = ['W1']
    species_lst = ['W1', 'P1']
    T_lst = [300.0, 350.0, 400.0, 450.0, 500.0, 550.0, 600.0,
             650.0, 700.0, 750.0, 800.0, 850.0, 900.0, 950.0, 1000.0,
             1050.0, 1100.0, 1150.0, 1200.0, 1250.0, 1300.0, 1350.0,
             1400.0, 1450.0, 1500.0, 1550.0, 1600.0, 1650.0, 1700.0,
             1750.0, 1800.0, 1850.0, 1900.0, 1950.0, 2000.0, 2050.0,
             2100.0, 2150.0, 2200.0, 2250.0, 2300.0, 2350.0, 2400.0,
             2450.0, 2500.0]
    P_lst = [0.01, 0.1, 0.316, 1.0, 3.16, 10.0, 31.6, 100.0]
    hoten = mess_io.reader.hoten.extract_hot_branching(HOT_OUT_SINGLE,
                                                     hotspecies_lst, species_lst, T_lst, P_lst)['W1']

    assert np.allclose(hoten[3.16][1200].iloc[0].values, np.array([1, 0]))
    assert np.allclose(hoten[3.16][1200].iloc[-1].values, np.array([0, 1]))
    assert np.allclose(hoten[3.16][1200].loc[21.1231].values, np.array([6.74809675e-05, 9.99932519e-01]))

    hotspecies_lst = ['CH3CH2CH2', 'CH3CHCH3']
    species_lst = ['CH3CH2CH2', 'CH3CHCH3', 'P1', 'P2']
    T_lst = [400.0, 450.0, 500.0, 550.0, 600.0, 650.0, 700.0, 750.0,
             800.0, 850.0, 900.0, 950.0, 1000.0, 1050.0, 1100.0, 1150.0,
             1200.0, 1250.0, 1300.0, 1350.0, 1400.0, 1450.0, 1500.0,
             1550.0, 1600.0, 1650.0, 1700.0, 1750.0, 1800.0, 1850.0,
             1900.0, 1950.0, 2000.0]
    P_lst = [0.1, 1.0, 10.0, 100.0]
    hoten = mess_io.reader.hoten.extract_hot_branching(HOT_OUT_DOUBLE,
                                                     hotspecies_lst, species_lst, T_lst, P_lst)['CH3CHCH3']
    assert np.allclose(hoten[100][1800].iloc[0].values, np.array([0.0, 1.000000, 0.000000, 0.000000]))
    assert np.allclose(hoten[100][1800].iloc[-1].values, np.array([0.0, 0.000117, 0.046422, 0.953461]), atol=1e-5)
    assert np.allclose(hoten[100][1800].loc[149.174].values, np.array([0.0, 0.000224, 0.045513, 0.954262]), atol=1e-5)

if __name__ == '__main__':
    test_get_hot_names()
    test_extract_hot_branching()
