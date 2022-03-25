"""
  Convert the PAC99 polynomial
"""


def pac2ckin_poly(name, atom_dict, pac99_poly_str):
    """ Convert a NASA polynomial from its format given in PAC99 input
        to the format used in ChemKin mechanisms.

        :param name: species name
        :type name: str
        :param atom_dict: number of each atom in the species
        :type atom_dict: dict[str: int]
        :param pac99_poly_str: PAC99-format NASA polynomial
        :type pac99_poly_str: str
        :rtype: str
    """

    # Parse the lines of the pac string containing the desired coefficients
    las = [0.0 for i in range(7)]
    has = [0.0 for i in range(7)]
    las[0:5] = _parse_line16(pac99_poly_str.splitlines()[6][0:80])
    las[5:7] = _parse_line16(pac99_poly_str.splitlines()[7][48:80])
    has[0:5] = _parse_line16(pac99_poly_str.splitlines()[9][0:80])
    has[5:7] = _parse_line16(pac99_poly_str.splitlines()[10][48:80])

    if 'H' in atom_dict:
        num_h = atom_dict['H']
    else:
        num_h = 0
    if 'C' in atom_dict:
        num_c = atom_dict['C']
    else:
        num_c = 0
    if 'N' in atom_dict:
        num_n = atom_dict['N']
    else:
        num_n = 0
    if 'O' in atom_dict:
        num_o = atom_dict['O']
    else:
        num_o = 0
    if 'F' in atom_dict:
        num_f = atom_dict['F']
    else:
        num_f = 0

    # Build a string for the NASA polynomial in ChemKin format
    line1 = "%s        H %3dC %3dO %3dN %3dG%10.1F%10.1F%8.1FF %3d 1\n" % (
        name.ljust(16)[0:16],
        num_h, num_c, num_o, num_n, 
        200.0, 3000.0, 1000.0, num_f)
    line2 = "% 15.8E% 15.8E% 15.8E% 15.8E% 15.8E    2\n" % (
        has[0], has[1], has[2], has[3], has[4])
    line3 = "% 15.8E% 15.8E% 15.8E% 15.8E% 15.8E    3\n" % (
        has[5], has[6], las[0], las[1], las[2])
    line4 = "% 15.8E% 15.8E% 15.8E% 15.8E                   4\n" % (
        las[3], las[4], las[5], las[6])

    return line1 + line2 + line3 + line4


def _parse_line16(string):
    """ Parse the values of a string containing exponental numbers
        of length 16 chars, where the numbers may not have a space in between.

        :param string: string containing the exponential numbers
        :type string: str
        :rtype: list(float)
    """

    assert len(string) % 16 == 0, 'Given string should have 16n chararacters'

    # Replace the exponent D with E
    string2 = string.replace('D', 'E')

    # Build a list of values from the string
    nchunks = len(string2) // 16
    vals = [0.0 for i in range(nchunks)]
    for i in range(nchunks):
        vals[i] = float(string2[i*16: (i+1)*16])

    return vals
