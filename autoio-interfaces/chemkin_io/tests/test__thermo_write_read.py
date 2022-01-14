""" test chemkin_io.writer.mechanism.thermo_block
"""

from chemkin_io.writer.mechanism import thermo_block as writer
from chemkin_io.parser.thermo import create_spc_nasa7_dct as parser


# Define a sample spc_nasa7_dct
SPC_NASA7_DCT = {
    'O2': ['RUS 89', 'O   2               ', 'G', [200.0, 6000.0, 1000.0],
           ([2.54363697E+00, -2.73162486E-05, -4.19029520E-09, 4.95481845E-12,
             -4.79553694E-16, 2.92260120E+04, 4.92229457E+00],
            [3.16826710E+00, -3.27931884E-03, 6.64306396E-06, -6.12806624E-09,
             2.11265971E-12, 2.91222592E+04, 2.05193346E+00])
           ],
    'N2O': ['L 7/88', 'N   1O   1          ', 'G', [200.0, 6000.0, 1000.0],
            ([0.48230729E+01, 0.26270251E-02, -0.95850872E-06, 0.16000712E-09,
              -0.97752302E-14, 0.80734047E+04, -0.22017208E+01],
             [0.22571502E+01, 0.11304728E-01, -0.13671319E-04, 0.96819803E-08,
              -0.29307182E-11, 0.87417746E+04, 0.10757992E+02])
            ]
}


def test__thermo_read_write():
    """ Tests the chemkin_io parsing and writing for thermo
    """

    # Run the writer
    thermo_str = writer(SPC_NASA7_DCT)

    # Run the parser
    new_spc_nasa7_dct = parser(thermo_str)

    # Run the writer again on the new spc_nasa7_dct
    new_thermo_str = writer(new_spc_nasa7_dct)

    assert thermo_str == new_thermo_str
