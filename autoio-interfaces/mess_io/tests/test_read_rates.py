""" test mess_io.reader.rates
"""

import os
import numpy
import automol.util.dict_
from ioformat import pathtools
import mess_io.reader


PATH = os.path.dirname(os.path.realpath(__file__))
INP_PATH = os.path.join(PATH, 'data', 'inp')
OUT_PATH = os.path.join(PATH, 'data', 'out')

KTP_INP_STR = pathtools.read_file(INP_PATH, 'example.inp')
KTP_OUT_STR = pathtools.read_file(OUT_PATH, 'rate.out')
KTP_OUT_BAR_STR = pathtools.read_file(OUT_PATH, 'rate.out_bar')
KTP_OUT_TORR_STR = pathtools.read_file(OUT_PATH, 'rate.out_torr')
KE_OUT_STR = pathtools.read_file(OUT_PATH, 'ke.out')

# Set the REACTANT and PRODUCT
REACTANT = 'W1'
PRODUCT = 'W3'

# Define things for testing the filtering function
KTP_DCT1 = { 
    # Weird negative, non-physical k(T)s and undefined k(T)s
    0.01: ((200., 400., 600., 800., 1000., 1200.),
           (-1.0e3, 2.0e3, -3.0e3, None, None, None)),
    # Just negative k(T) at first position and undefined k(T)s
    0.1: ((200., 400., 600., 800., 1000., 1200.),
          (-1.0e4, 2.0e4, 3.0e4, None, None, None)),
    # Positive k(T) and undefined k(T)s
    1.0: ((200., 400., 600., 800., 1000., 1200.),
          (1.0e5, 2.0e5, 3.0e5, None, None, None)),
    # Positive k(T)s
    10.0: ((200., 400., 600., 800., 1000., 1200.),
           (1.0e6, 2.0e6, 3.0e6, 4.0e6, 5.0e6, 6.0e6)),
    # Positive k(T)s
    100.0: ((200., 400., 600., 800., 1000., 1200.),
            (1.0e7, 2.0e7, 3.0e7, 4.0e7, 5.0e7, 6.0e7))}
# Handles very small rate constants that may not wish to be modelled
KTP_DCT2 = { 
    0.01: ((200., 400., 600., 800., 1000., 1200.),
           (1.0e-25, 5.0e-25, 9.0e-24, 1.5e-23, 6.5e-22, 1.0e-20)),
    1.0: ((200., 400., 600., 800., 1000., 1200.),
          (8.0e-25, 3.0e-24, 2.0e-22, 6.0e-20, 2.5e-18, 5.0e-16)),}
# Dict that should come back empty after filter
KTP_DCT3 = { 
    0.01: ((200., 400., 600., 800., 1000., 1200.),
           (-1.0e-25, -5.0e-25, -9.0e-24, -1.5e-23, -6.5e-22, -1.0e-20)),
    1.0: ((200., 400., 600., 800., 1000., 1200.),
          (-8.0e-25, -3.0e-24, -2.0e-22, -6.0e-20, -2.5e-18, -5.0e-16)),}


def test__ktp_dct():
    """ test mess_io.reader.rates.ktp_dct
    """

    ref_ktp_dct = {
        'high': (
            (500.0, 650.0, 800.0, 950.0,
             1100.0, 1250.0, 1400.0, 1550.0,
             1700.0, 1850.0, 2000.0),
            (0.00123704, 2.07339, 295.551, 9910.26,
             135400.0, 1024120.0, 5143860.0, 19234000.0,
             57729800.0, 146453000.0, 325721000.0)),
        0.1: (
            (500.0, 650.0, 800.0, 950.0,
             1100.0, 1250.0, 1400.0, 1550.0,
             1700.0, 1850.0, 2000.0),
            (0.000536676, 0.127095, 3.43383, 25.6351,
             84.5986, 327.575, 1392.42, 4112.93,
             9141.18, 14058.6, 19558.3)),
        1.0: (
            (500.0, 650.0, 800.0, 950.0,
             1100.0, 1250.0, 1400.0, 1550.0,
             1700.0, 1850.0, 2000.0),
            (0.000892079, 0.395065, 11.9261, 114.45,
             473.081, 1469.97, 6064.93, 21266.3,
             62019.4, 138742.0, 279227.0)),
        10.0: (
            (500.0, 650.0, 800.0, 950.0,
             1100.0, 1250.0, 1400.0, 1550.0,
             1700.0, 1850.0, 2000.0),
            (0.00114645, 1.04871, 46.0872, 471.742,
             2410.59, 7623.67, 23580.8, 76336.0,
             235593.0, 587439.0, 1322680.0)),
        100.0: (
            (500.0, 650.0, 800.0, 950.0,
             1100.0, 1250.0, 1400.0, 1550.0,
             1700.0, 1850.0, 2000.0),
            (0.00122221, 1.76651, 145.044, 1977.75,
             10753.1, 37345.5, 105181.0, 268989.0,
             728117.0, 1802380.0, 4143210.0))
    }

    ktp_dct = mess_io.reader.rates.ktp_dct(
        KTP_OUT_STR, REACTANT, PRODUCT)

    assert set(ktp_dct.keys()) == set(ref_ktp_dct.keys())
    for pressure, tk_arr in ktp_dct.items():
        assert numpy.allclose(tk_arr, ref_ktp_dct[pressure])

    # Read files that have units that are in bar, torr instead of atm
    # vals should be same as above; just changed units in output string
    # for testing
    ktp_dct_bar = mess_io.reader.rates.ktp_dct(
        KTP_OUT_BAR_STR, REACTANT, PRODUCT)
    ktp_dct_torr = mess_io.reader.rates.ktp_dct(
        KTP_OUT_TORR_STR, REACTANT, PRODUCT)

    tkbarr = automol.util.dict_.value_in_floatkey_dct(
        ktp_dct_bar, 0.98692, tol=0.01)
    tktorr = automol.util.dict_.value_in_floatkey_dct(
        ktp_dct_torr, 0.0013157894736842107, tol=0.00001)

    assert numpy.allclose(ref_ktp_dct[1.0], tkbarr)
    assert numpy.allclose(ref_ktp_dct[1.0], tktorr)


def test__ke_dct():
    """ test mess_io.reader.rates.ke_dct
    """

    ref_ke_dct = {
        0.0: 0.0, 0.2: 7.76e-15, 0.4: 1.1e-13,
        0.6: 4.49e-13, 0.8: 1.54e-12,
        1.0: 4.53e-12, 1.2: 1.25e-11, 1.4: 3.15e-11,
        1.6: 7.32e-11, 1.8: 1.67e-10, 2.0: 3.66e-10,
        2.2: 7.76e-10, 2.4: 1.59e-09, 2.6: 3.22e-09,
        2.8: 6.36e-09, 3.0: 1.23e-08, 3.2: 2.36e-08,
        3.4: 4.44e-08, 3.6: 8.25e-08, 3.8: 1.51e-07,
        4.0: 2.75e-07, 4.2: 4.93e-07, 4.4: 8.77e-07,
        4.6: 1.55e-06, 4.8: 2.7e-06, 5.0: 4.68e-06,
        5.2: 8.06e-06, 5.4: 1.38e-05, 5.6: 2.34e-05,
        5.8: 3.94e-05, 6.0: 6.6e-05, 6.2: 0.000109997,
        6.4: 0.000182243, 6.6: 0.000300318, 6.8: 0.000492433,
        7.0: 0.000803405, 7.2: 0.0013047, 7.4: 0.00210909,
        7.6: 0.00339475, 7.8: 0.00544102, 8.0: 0.00868544,
        8.2: 0.0138106, 8.4: 0.0218763, 8.6: 0.0345266,
        8.8: 0.0542978, 9.0: 0.085099, 9.2: 0.132925,
        9.4: 0.206959, 9.6: 0.321211, 9.8: 0.497013,
        10.0: 0.766741, 10.2: 1.17942, 10.4: 1.80909,
        10.6: 2.76729, 10.8: 4.22164, 11.0: 6.42348,
        11.2: 9.74868, 11.4: 14.7583, 11.6: 22.2877,
        11.8: 33.5782, 12.0: 50.4699}

    ke_dct = mess_io.reader.rates.ke_dct(
        KE_OUT_STR, REACTANT, PRODUCT)

    assert set(ke_dct.keys()) == set(ref_ke_dct.keys())
    for ene, ratek in ke_dct.items():
        assert numpy.isclose(ratek, ref_ke_dct[ene])


def test__tp():
    """ test mess_io.reader.rates.pressures
        test mess_io.reader.rates.temperatures
    """

    ref_inp_temps = (600.0, 800.0, 1000.0, 1200.0, 1400.0, 1600.0, 1800.0,
                     2000.0, 2200.0, 2400.0, 2600.0, 2800.0, 3000.0)
    ref_out_temps = (500.0, 650.0, 800.0, 950.0, 1100.0, 1250.0, 1400.0,
                     1550.0, 1700.0, 1850.0, 2000.0)
    ref_inp_press = (0.01, 0.1, 1.0, 10.0, 100.0, 'high')
    ref_out_press = (0.1, 1.0, 10.0, 100.0, 'high')

    inp_temps, inp_tunit = mess_io.reader.rates.temperatures(
        KTP_INP_STR, mess_file='inp')
    out_temps, out_tunit = mess_io.reader.rates.temperatures(
        KTP_OUT_STR, mess_file='out')

    inp_press, inp_punit = mess_io.reader.rates.pressures(
        KTP_INP_STR, mess_file='inp')
    out_press, out_punit = mess_io.reader.rates.pressures(
        KTP_OUT_STR, mess_file='out')

    assert numpy.allclose(ref_inp_temps, inp_temps)
    assert numpy.allclose(ref_out_temps, out_temps)
    assert inp_tunit == out_tunit == 'K'

    assert ref_inp_press[-1] == inp_press[-1]
    assert ref_out_press[-1] == out_press[-1]
    assert numpy.allclose(ref_inp_press[:-1], inp_press[:-1])
    assert numpy.allclose(ref_out_press[:-1], out_press[:-1])
    assert inp_punit == out_punit == 'atm'


def test__rxns_labels():
    """ test mess_io.reader.rates.reactions
    """

    ref_rxns1 = (
        ('W1', 'W2'), ('W1', 'W3'), ('W1', 'W4'), ('W1', 'P1'), ('W1', 'P2'),
        ('W1', 'P3'), ('W1', 'P4'), ('W1', 'P5'), ('W1', 'P6'), ('W2', 'W1'),
        ('W2', 'W3'), ('W2', 'W4'), ('W2', 'P1'), ('W2', 'P2'), ('W2', 'P3'),
        ('W2', 'P4'), ('W2', 'P5'), ('W2', 'P6'), ('W3', 'W1'), ('W3', 'W2'),
        ('W3', 'W4'), ('W3', 'P1'), ('W3', 'P2'), ('W3', 'P3'), ('W3', 'P4'),
        ('W3', 'P5'), ('W3', 'P6'), ('W4', 'W1'), ('W4', 'W2'), ('W4', 'W3'),
        ('W4', 'P1'), ('W4', 'P2'), ('W4', 'P3'), ('W4', 'P4'), ('W4', 'P5'),
        ('W4', 'P6'), ('P1', 'W1'), ('P1', 'W2'), ('P1', 'W3'), ('P1', 'W4'),
        ('P1', 'P2'), ('P1', 'P3'), ('P1', 'P4'), ('P1', 'P5'), ('P1', 'P6'),
        ('P2', 'W1'), ('P2', 'W2'), ('P2', 'W3'), ('P2', 'W4'), ('P2', 'P1'),
        ('P2', 'P3'), ('P2', 'P4'), ('P2', 'P5'), ('P2', 'P6'), ('P3', 'W1'),
        ('P3', 'W2'), ('P3', 'W3'), ('P3', 'W4'), ('P3', 'P1'), ('P3', 'P2'),
        ('P3', 'P4'), ('P3', 'P5'), ('P3', 'P6'), ('P4', 'W1'), ('P4', 'W2'),
        ('P4', 'W3'), ('P4', 'W4'), ('P4', 'P1'), ('P4', 'P2'), ('P4', 'P3'),
        ('P4', 'P5'), ('P4', 'P6'), ('P5', 'W1'), ('P5', 'W2'), ('P5', 'W3'),
        ('P5', 'W4'), ('P5', 'P1'), ('P5', 'P2'), ('P5', 'P3'), ('P5', 'P4'),
        ('P5', 'P6'), ('P6', 'W1'), ('P6', 'W2'), ('P6', 'W3'), ('P6', 'W4'),
        ('P6', 'P1'), ('P6', 'P2'), ('P6', 'P3'), ('P6', 'P4'), ('P6', 'P5'))
    ref_rxns2 = (
        ('W1', 'W2'), ('W1', 'W3'), ('W1', 'W4'), ('W1', 'P1'), ('W1', 'P2'),
        ('W1', 'P3'), ('W1', 'P4'), ('W1', 'P5'), ('W1', 'P6'), ('W2', 'W3'),
        ('W2', 'W4'), ('W2', 'P1'), ('W2', 'P2'), ('W2', 'P3'), ('W2', 'P4'),
        ('W2', 'P5'), ('W2', 'P6'), ('W3', 'W4'), ('W3', 'P1'), ('W3', 'P2'),
        ('W3', 'P3'), ('W3', 'P4'), ('W3', 'P5'), ('W3', 'P6'), ('W4', 'P1'),
        ('W4', 'P2'), ('W4', 'P3'), ('W4', 'P4'), ('W4', 'P5'), ('W4', 'P6'),
        ('P1', 'P2'), ('P1', 'P3'), ('P1', 'P4'), ('P1', 'P5'), ('P1', 'P6'),
        ('P2', 'P3'), ('P2', 'P4'), ('P2', 'P5'), ('P2', 'P6'), ('P3', 'P4'),
        ('P3', 'P5'), ('P3', 'P6'), ('P4', 'P5'), ('P4', 'P6'), ('P5', 'P6'))

    assert ref_rxns1 == mess_io.reader.rates.reactions(
        KTP_OUT_STR, read_rev=True)
    assert ref_rxns2 == mess_io.reader.rates.reactions(
        KTP_OUT_STR, read_rev=False)


def test__filter_ktp():
    """ test mess_io.reader.rates.filter_ktp_dct
    """

    # Check the temperatures and undefined (these are wrong)
    ref_dct1 = {
        0.1: (numpy.array([400., 600.]),
              numpy.array([2.0e4, 3.0e4])),
        1.0: (numpy.array([200., 400., 600.]),
              numpy.array([1.0e5, 2.0e5, 3.0e5])),
        10.0: (numpy.array([200.,  400.,  600.,  800., 1000., 1200.]),
               numpy.array([1.0e6, 2.0e6, 3.0e6, 4.0e6, 5.0e6, 6.0e6])),
        100.0: (numpy.array([200.,  400.,  600.,  800., 1000., 1200.]),
                numpy.array([1.0e7, 2.0e7, 3.0e7, 4.0e7, 5.0e7, 6.0e7]))}
    ref_dct2 = {
        0.01: (numpy.array([400.]),
               numpy.array([2.0e3])),
        0.1: (numpy.array([400., 600.]),
              numpy.array([2.0e4, 3.0e4])),
        1.0: (numpy.array([400., 600.]),
              numpy.array([2.0e5, 3.0e5])),
        10.0: (numpy.array([400., 600., 800.]),
               numpy.array([2.0e6, 3.0e6, 4.0e6])),
        100.0: (numpy.array([400., 600., 800.]),
                numpy.array([2.0e7, 3.0e7, 4.0e7]))}

    filt_ktp_dct_temp1 = mess_io.reader.rates.filter_ktp_dct(
        KTP_DCT1, bimol=False, tmin=None, tmax=None)

    filt_ktp_dct_temp2 = mess_io.reader.rates.filter_ktp_dct(
        KTP_DCT1, bimol=False, tmin=400.0, tmax=800.0)

    for key1, key2 in zip(filt_ktp_dct_temp1.keys(), ref_dct1.keys()):
        assert numpy.isclose(key1, key2)
        assert numpy.allclose(filt_ktp_dct_temp1[key1][0], ref_dct1[key2][0])
        assert numpy.allclose(filt_ktp_dct_temp1[key1][1], ref_dct1[key2][1])

    for key1, key2 in zip(filt_ktp_dct_temp2.keys(), ref_dct2.keys()):
        assert numpy.isclose(key1, key2)
        assert numpy.allclose(filt_ktp_dct_temp2[key1][0], ref_dct2[key2][0])
        assert numpy.allclose(filt_ktp_dct_temp2[key1][1], ref_dct2[key2][1])

    # Check for small rate constants (based on bimol flag)
    ref_dct3 = {
        0.01: (
            numpy.array([200.,  400.,  600.,  800., 1000., 1200.]),
            numpy.array([1.0e-25, 5.0e-25, 9.0e-24, 1.5e-23, 6.5e-22, 1.0e-20])),
        1.0: (
            numpy.array([200.,  400.,  600.,  800., 1000., 1200.]),
            numpy.array([8.0e-25, 3.0e-24, 2.0e-22, 6.0e-20, 2.5e-18, 5.0e-16]))}
    ref_dct4 = {
        0.01: (numpy.array([600.,  800., 1000., 1200.]),
               numpy.array([9.0e-24, 1.5e-23, 6.5e-22, 1.0e-20])),
        1.0: (numpy.array([400.,  600.,  800., 1000., 1200.]),
              numpy.array([3.0e-24, 2.0e-22, 6.0e-20, 2.5e-18, 5.0e-16]))}

    filt_ktp_dct_smallk1 = mess_io.reader.rates.filter_ktp_dct(
        KTP_DCT2, bimol=False, tmin=None, tmax=None)
    filt_ktp_dct_smallk2 = mess_io.reader.rates.filter_ktp_dct(
        KTP_DCT2, bimol=True, tmin=None, tmax=None)

    for key1, key2 in zip(filt_ktp_dct_smallk1.keys(), ref_dct3.keys()):
        assert numpy.isclose(key1, key2)
        assert numpy.allclose(filt_ktp_dct_smallk1[key1][0], ref_dct3[key2][0])
        assert numpy.allclose(filt_ktp_dct_smallk1[key1][1], ref_dct3[key2][1])

    for key1, key2 in zip(filt_ktp_dct_smallk2.keys(), ref_dct4.keys()):
        assert numpy.isclose(key1, key2)
        assert numpy.allclose(filt_ktp_dct_smallk2[key1][0], ref_dct4[key2][0])
        assert numpy.allclose(filt_ktp_dct_smallk2[key1][1], ref_dct4[key2][1])

    # Dict that should come back empty
    filt_ktp_dct_emptyret = mess_io.reader.rates.filter_ktp_dct(
        KTP_DCT3, bimol=False, tmin=None, tmax=None)

    assert not filt_ktp_dct_emptyret

    # Test the removal of pressures
    ref_dct5 = {
        0.1: (numpy.array([400., 600.]),
              numpy.array([2.0e4, 3.0e4])),
        1.0: (numpy.array([200., 400., 600.]),
              numpy.array([1.0e5, 2.0e5, 3.0e5])),
        10.0: (numpy.array([200.,  400.,  600.,  800., 1000., 1200.]),
               numpy.array([1.0e6, 2.0e6, 3.0e6, 4.0e6, 5.0e6, 6.0e6]))}
    filt_ktp_dct_pressure = mess_io.reader.rates.filter_ktp_dct(
        KTP_DCT1, bimol=False, tmin=None, tmax=None, pmin=0.1, pmax=10)
    for key1, key2 in zip(filt_ktp_dct_pressure.keys(), ref_dct5.keys()):
        assert numpy.isclose(key1, key2)
        assert numpy.allclose(filt_ktp_dct_pressure[key1][0], ref_dct5[key2][0])
        assert numpy.allclose(filt_ktp_dct_pressure[key1][1], ref_dct5[key2][1])


if __name__ == '__main__':
    test__ktp_dct()
    test__ke_dct()
    test__tp()
    test__rxns_labels()
    test__filter_ktp()
