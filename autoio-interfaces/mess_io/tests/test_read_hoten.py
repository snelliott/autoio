""" test mess_io.reader.hoten
"""

import os
import numpy as np
import ioformat
import mess_io


PATH = os.path.dirname(os.path.realpath(__file__))
IPATH = os.path.join(PATH, 'data', 'inp')
OPATH = os.path.join(PATH, 'data', 'out')
HOT_INP_SGL = ioformat.pathtools.read_file(IPATH, 'me_ktp_hoten_ch2o_oh.inp')
HOT_LOG_SGL = ioformat.pathtools.read_file(OPATH, 'me_ktp_hoten_ch2o_oh.logf')
HOT_INP_DBL = ioformat.pathtools.read_file(IPATH, 'me_ktp_hoten_c3h7.inp')
HOT_LOG_DBL = ioformat.pathtools.read_file(OPATH, 'me_ktp_hoten_c3h7.logf')


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
        HOT_LOG_SGL, hotspecies_en, species_lst, temp_lst, pressure_lst)
    hoten = hoten_branch_dct['W1']

    assert np.allclose(hoten[3.16][1200].iloc[0].values,
                       np.array([0.99892069, 0.00107931]))
    assert np.allclose(hoten[3.16][1200].iloc[-1].values,
                       np.array([0, 1]))
    assert np.allclose(hoten[3.16][1200].loc[21.1231].values,
                       np.array([6.74809675e-05, 9.99932519e-01]))

    hotspecies_en = {'CH3CH2CH2': 3.19, 'CH3CHCH3': 0.0}
    species_lst = ('CH3CH2CH2', 'CH3CHCH3', 'P1', 'P2')
    temp_lst = (400.0, 450.0, 500.0, 550.0, 600.0, 650.0, 700.0, 750.0,
                800.0, 850.0, 900.0, 950.0, 1000.0, 1050.0, 1100.0, 1150.0,
                1200.0, 1250.0, 1300.0, 1350.0, 1400.0, 1450.0, 1500.0,
                1550.0, 1600.0, 1650.0, 1700.0, 1750.0, 1800.0, 1850.0,
                1900.0, 1950.0, 2000.0)
    pressure_lst = (0.1, 1.0, 10.0, 100.0)
    hoten_branch_dct = mess_io.reader.hoten.extract_hot_branching(
        HOT_LOG_DBL, hotspecies_en, species_lst, temp_lst, pressure_lst)
    hoten = hoten_branch_dct['CH3CHCH3']

    assert np.allclose(hoten[100][1800].iloc[0].values,
                       np.array([3.05556914e-04, 9.97009793e-01, 0.00000000e+00, 2.68464971e-03]))
    assert np.allclose(hoten[100][1800].iloc[-1].values,
                       np.array([0.00000000e+00, 1.17216104e-04, 4.64219410e-02, 9.53460843e-01]),
                       atol=1e-5)
    assert np.allclose(hoten[100][1800].loc[149.174].values,
                       np.array([0.00000000e+00, 2.24250549e-04, 4.55134115e-02, 9.54262338e-01]),
                       atol=1e-5)

if __name__ == '__main__':
    test_get_hot_species()
    test_extract_hot_branching()