""" test autorun.nst
"""

import os
# import tempfile
import automol
from ioformat import pathtools
import autorun


PATH = os.path.dirname(os.path.realpath(__file__))
DAT_PATH = os.path.join(PATH, 'data')


# Common structural info
CHARGE = 0
MULTS = (1, 3)

# NNO - MOLPRO
NNO_GEO = automol.geom.from_string(
    pathtools.read_file(DAT_PATH, 'nno_init.xyz'))
NNO_PROG = 'molpro2015'
NNO_METHOD = 'caspt2'
NNO_BASIS = 'aug-cc-pvdz'
NNO_ORB_TYPE = 'RR'
NNO_INI_KWARGS = {
    'casscf_options': ('closed=5', 'wf,22,1,0,0'),
    'core_options': ('shift=0.20',)
}

# C2H4O - G09
C2H4O_GEO = automol.geom.from_string(
    pathtools.read_file(DAT_PATH, 'c2h4o_init.xyz'))

C2H4O_PROG = 'gaussian09'
C2H4O_METHOD = 'm062x'
C2H4O_BASIS = 'cc-pvdz'
C2H4O_ORB_TYPE = 'RU'
C2H4O_INI_KWARGS = {
    'gen_lines': {1: ('# int=ultrafine',)},
    'machine_options': ('%NProcShared=10',)
}


def test__molpro_flux():
    """ test autorun.nst
    """

    run_dir = 'tmp-mol'
    # with tempfile.TemporaryDirectory(dir=PATH) as run_dir:
    isc = autorun.nst.isc_flux(
        run_dir,
        NNO_PROG, NNO_GEO, CHARGE, MULTS,
        NNO_METHOD, NNO_BASIS, NNO_ORB_TYPE, NNO_INI_KWARGS)
    print(isc)


def test__gaussian_flux():
    """ test autorun.nst
    """

    run_dir = 'tmp-g09'
    # with tempfile.TemporaryDirectory(dir=PATH) as run_dir:
    isc = autorun.nst.isc_flux(
        run_dir,
        C2H4O_PROG, C2H4O_GEO, CHARGE, MULTS,
        C2H4O_METHOD, C2H4O_BASIS, C2H4O_ORB_TYPE, C2H4O_INI_KWARGS)
    print(isc)


if __name__ == '__main__':
    # test__molpro_flux()
    test__gaussian_flux()
