""" test mess_io.reader.ped
"""

import os
import numpy as np
import ioformat
import mess_io


PATH = os.path.dirname(os.path.realpath(__file__))
INP_PATH_SGL = os.path.join(PATH, 'data', 'prompt', 'CH2O_H')
INP_PATH_DBL = os.path.join(PATH, 'data', 'prompt', 'C3H8_H')
PED_INP_SGL = ioformat.pathtools.read_file(INP_PATH_SGL, 'me_ktp_ped.inp')
PED_OUT_SGL = ioformat.pathtools.read_file(INP_PATH_SGL, 'thf_pyro_ped.out')
PED_INP_DBL = ioformat.pathtools.read_file(INP_PATH_DBL, 'me_ktp_ped.inp')
PED_OUT_DBL = ioformat.pathtools.read_file(INP_PATH_DBL, 'thf_pyro_ped.out')


def test_ped_names():
    """ test mess_io.reader.ped.ped_names
    """

    pedspecies, pedoutput = mess_io.reader.ped.ped_names(PED_INP_SGL)
    assert pedspecies == (('RH', 'R',),)
    assert pedoutput == 'thf_pyro_ped.out'

    pedspecies, pedoutput = mess_io.reader.ped.ped_names(PED_INP_DBL)
    assert pedspecies == (('RH', 'NC3H7'), ('RH', 'IC3H7'))
    assert pedoutput == 'thf_pyro_ped.out'


def test_ped_get_ped():
    """ test mess_io.reader.ped.ped_names
    """

    # PES1
    pedspecies1 = (('RH', 'R'),)
    energy_dct1 = {'W0': -0.8, 'RH': 0.0, 'R': -16.7,
                   'B0': 0.0, 'B1': 5.2}
    ped_dct1 = mess_io.reader.ped.get_ped(
        PED_OUT_SGL, pedspecies1, energy_dct1)

    assert np.isclose(
        ped_dct1['RH->R'][1][500].index[1], 16.814,
        atol=1e-3, rtol=1e-3)

    # PES2
    pedspecies2 = (('RH', 'NC3H7'), ('RH', 'IC3H7'))
    energy_dct2 = {'W0': 0.0, 'RH': 0.0, 'NC3H7': -3.53, 'IC3H7': -6.58,
                   'B0': 0.0, 'B1': 10.19, 'B2': 7.67}
    ped_dct2 = mess_io.reader.ped.get_ped(
        PED_OUT_DBL, pedspecies2, energy_dct2)

    print(ped_dct2['RH->NC3H7'][1][500].index[1])

    assert np.isclose(
        ped_dct2['RH->NC3H7'][1][500].index[1], 3.730,
        atol=1e-3, rtol=1e-3)
    assert np.isclose(
        ped_dct2['RH->IC3H7'][1][500].index[1], 6.780,
        atol=1e-3, rtol=1e-3)


if __name__ == '__main__':
    # test_ped_names()
    test_ped_get_ped()
