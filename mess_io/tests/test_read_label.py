""" test mess_io.reader.pes
"""

import os
from ioformat import pathtools
from mess_io.reader import relabel
from mess_io.reader import name_label_dct



PATH = os.path.dirname(os.path.realpath(__file__))
OPATH = os.path.join(PATH, 'data', 'out')
HOT_LOG_DBL = pathtools.read_file(OPATH, 'me_ktp_hoten_c3h7.logf')

LBL_DCT = {'W1': 'CH3CH2CH2', 'W2': 'CH3CHCH3', 'P1': 'C2H4+CH3', 'P2': 'CH3CHCH2+H'}

def test_name_label_dct():
    """ test mess_io._label.name_label_dct
    """
    lbl_dct = name_label_dct(HOT_LOG_DBL)
    assert lbl_dct == LBL_DCT


if __name__ == '__main__':
    test_name_label_dct()
    # test_relabel() #todo
