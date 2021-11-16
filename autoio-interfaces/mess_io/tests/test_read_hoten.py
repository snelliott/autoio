""" test mess_io.reader.hoten
"""

import os
import numpy as np
import ioformat
import mess_io


PATH = os.path.dirname(os.path.realpath(__file__))
INP_PATH_SGL = os.path.join(PATH, 'data', 'prompt', 'CH2O_OH_ORIG')
INP_PATH_DBL = os.path.join(PATH, 'data', 'prompt', 'C3H7')
HOT_INP_SGL = ioformat.pathtools.read_file(INP_PATH_SGL, 'me_ktp_hoten.inp')
HOT_LOG_SGL = ioformat.pathtools.read_file(INP_PATH_SGL, 'me_ktp_hoten.log')
HOT_INP_DBL = ioformat.pathtools.read_file(INP_PATH_DBL, 'me_ktp_hoten.inp')
HOT_LOG_DBL = ioformat.pathtools.read_file(INP_PATH_DBL, 'me_ktp_hoten.log')


def test_get_hot_names():
    """ test mess_io.read.get_hot_names
    """
    assert mess_io.reader.hoten.get_hot_names(HOT_INP_SGL) == (
        'W1',)
    assert mess_io.reader.hoten.get_hot_names(HOT_INP_DBL) == (
        'CH3CH2CH2', 'CH3CHCH3')


def test_extract_hot_branching():
    """ test mess_io.read.extract_hot_branching
    """
    hotspecies_lst = ('W1',)
    species_lst = ('W1', 'P1')
    temp_lst = (300.0, 350.0, 400.0, 450.0, 500.0, 550.0, 600.0,
                650.0, 700.0, 750.0, 800.0, 850.0, 900.0, 950.0, 1000.0,
                1050.0, 1100.0, 1150.0, 1200.0, 1250.0, 1300.0, 1350.0,
                1400.0, 1450.0, 1500.0, 1550.0, 1600.0, 1650.0, 1700.0,
                1750.0, 1800.0, 1850.0, 1900.0, 1950.0, 2000.0, 2050.0,
                2100.0, 2150.0, 2200.0, 2250.0, 2300.0, 2350.0, 2400.0,
                2450.0, 2500.0)
    pressure_lst = (0.01, 0.1, 0.316, 1.0, 3.16, 10.0, 31.6, 100.0)
    hoten_branch_dct = mess_io.reader.hoten.extract_hot_branching(
        HOT_LOG_SGL, hotspecies_lst, species_lst, temp_lst, pressure_lst)
    hoten = hoten_branch_dct['W1']

    assert np.allclose(hoten[3.16][1200].iloc[0].values,
                       np.array([0, 1]))
    assert np.allclose(hoten[3.16][1200].iloc[-1].values,
                       np.array([1, 0]))
    assert np.allclose(hoten[3.16][1200].loc[21.1231].values,
                       np.array([6.74809675e-05, 9.99932519e-01]))

    hotspecies_lst = ('CH3CH2CH2', 'CH3CHCH3')
    species_lst = ('CH3CH2CH2', 'CH3CHCH3', 'P1', 'P2')
    temp_lst = (400.0, 450.0, 500.0, 550.0, 600.0, 650.0, 700.0, 750.0,
                800.0, 850.0, 900.0, 950.0, 1000.0, 1050.0, 1100.0, 1150.0,
                1200.0, 1250.0, 1300.0, 1350.0, 1400.0, 1450.0, 1500.0,
                1550.0, 1600.0, 1650.0, 1700.0, 1750.0, 1800.0, 1850.0,
                1900.0, 1950.0, 2000.0)
    pressure_lst = (0.1, 1.0, 10.0, 100.0)
    hoten_branch_dct = mess_io.reader.hoten.extract_hot_branching(
        HOT_LOG_DBL, hotspecies_lst, species_lst, temp_lst, pressure_lst)
    hoten = hoten_branch_dct['CH3CHCH3']

    print(hoten)
    assert np.allclose(hoten[100][1800].iloc[0].values,
                       np.array([0.0, 1.000000, 0.000000, 0.000000]))
    assert np.allclose(hoten[100][1800].iloc[-1].values,
                       np.array([0.0, 0.000117, 0.046422, 0.953461]),
                       atol=1e-5)
    assert np.allclose(hoten[100][1800].loc[149.174].values,
                       np.array([0.0, 0.000224, 0.045513, 0.954262]),
                       atol=1e-5)


if __name__ == '__main__':
    test_get_hot_names()
    test_extract_hot_branching()
