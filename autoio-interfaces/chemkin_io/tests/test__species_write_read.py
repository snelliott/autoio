""" tests chemkin_io.writer.mechanism.species_block
"""

from chemkin_io.writer.mechanism import species_block as writer
from chemkin_io.parser.species import names as parser


SPC_IDENT_DCT = {
    'O': {'smiles': 'smiles_1',
          'inchi': 'inchi_1',
          'charge': '',
          'mult': '',
          'sens': ''},
    'H': {'smiles': 'smiles_2',
          'inchi': 'inchi_2',
          'charge': '',
          'mult': '',
          'sens': ''}
}
SPC_NAMES_STRING_1 = (
    'SPECIES \n\nO     ! SMILES: smiles_1       ' +
    'InChi: inchi_1  \nH     ! SMILES: smiles_2       ' +
    'InChi: inchi_2  \n\nEND \n\n\n'
)
SPC_NAMES_STRING_2 = 'OH \nHO2 \nC3H8 \nN2O'
SPC_NAMES_TUPLE = ('OH', 'HO2', 'C3H8', 'N2O')


def test__write_spc_names():
    """ Tests the species names writing
    """
    spc_str = writer(SPC_IDENT_DCT)
    assert spc_str == SPC_NAMES_STRING_1


def test__read_spc_names():
    """ Tests the parsing of species names
    """
    spc_tuple = parser(SPC_NAMES_STRING_2)
    assert spc_tuple == SPC_NAMES_TUPLE
