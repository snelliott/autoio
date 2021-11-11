""" test mess_io.reader.ped
"""

import os
import numpy as np
import mess_io
from autofile.io_ import read_file

PATH = os.path.dirname(os.path.realpath(__file__))
INP_PATH_SINGLE = os.path.join(PATH, 'test_nB', 'CH2O_H')
PED_INP_SINGLE = read_file(os.path.join(INP_PATH_SINGLE, 'me_ktp_ped.inp'))
PED_OUT_SINGLE = read_file(os.path.join(INP_PATH_SINGLE, 'thf_pyro_ped.out'))
INP_PATH_DOUBLE = os.path.join(PATH, 'test_nB', 'C3H8_H')
PED_INP_DOUBLE = read_file(os.path.join(INP_PATH_DOUBLE, 'me_ktp_ped.inp'))
PED_OUT_DOUBLE = read_file(os.path.join(INP_PATH_DOUBLE, 'thf_pyro_ped.out'))

def test_ped_names():

    pedspecies, pedoutput = mess_io.reader.ped.ped_names(PED_INP_SINGLE)
    assert pedspecies == [['RH','R']]
    assert pedoutput == 'thf_pyro_ped.out'
    pedspecies, pedoutput = mess_io.reader.ped.ped_names(PED_INP_DOUBLE)
    assert pedspecies == [['RH', 'NC3H7'], ['RH', 'IC3H7']]
    assert pedoutput == 'thf_pyro_ped.out'

def test_ped_get_ped():

    # PES1
    pedspecies1 = [['RH', 'R']]
    energy_dct1 = {'W0': -0.8, 'RH': 0.0, 'R': -16.7, 'B0': 0.0, 'B1': 5.2}
    ped_dct1 = mess_io.reader.ped.get_ped(PED_OUT_SINGLE, pedspecies1, energy_dct1)
    assert np.isclose(ped_dct1['RH->R'][1][500].index[1], 16.814, atol=1e-3, rtol=1e-3)
    # PES2 
    pedspecies2 = [['RH', 'NC3H7'], ['RH', 'IC3H7']] 
    energy_dct2 = {'W0': 0.0, 'RH': 0.0, 'NC3H7': -3.53, 'IC3H7': -6.58, 'B0': 0.0, 'B1': 10.19, 'B2': 7.67}
    ped_dct2 = mess_io.reader.ped.get_ped(PED_OUT_DOUBLE, pedspecies2, energy_dct2)
    assert np.isclose(ped_dct2['RH->NC3H7'][1][500].index[1], 3.730, atol=1e-3, rtol=1e-3)
    assert np.isclose(ped_dct2['RH->IC3H7'][1][500].index[1], 6.780, atol=1e-3, rtol=1e-3)

if __name__ == '__main__':
    test_ped_names()
    test_ped_get_ped()

