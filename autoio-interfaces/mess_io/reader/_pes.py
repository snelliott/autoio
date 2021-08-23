"""
  Read a MESS input file and compile data for the PES inside
"""


def pes(input_string, read_fake=False):
    """ Read a MESS input file string and get info about PES

        :param input_string: string for a MESS (rates) input file
        :type input_string: str
        :param read_fake: value to include fake wells and barriers
        :type read_fake: bool
        :return energy_dct: dict[label: energy]
        :rtype: dict[label: energy]
        :return conn_lst
        :rtype: lst(str)
    """

    # Initialize energy and connection information
    energy_dct = {}
    conn_lst = tuple()
    pes_label_dct = {}

    input_lines = input_string.splitlines()
    for idx, line in enumerate(input_lines):

        if 'Well' in line:

            line_lst = line.strip().split()
            if line_lst[0].strip() == 'Well' and '!' not in line:
                # Get label
                label = line_lst[1]

                if ('F' not in label) or ('F' in label and read_fake):
                    # Get energy
                    for line2 in input_lines[idx:]:
                        if 'ZeroEnergy' in line2:
                            ene = float(line2.split()[-1])
                            break

                    # Add value to energy dct
                    energy_dct[label] = ene

                    # Add value to PES dct
                    # prior_line = input_lines[idx-1]
                    # line_lst2 = prior_line.split('!')
                    # spc = line_lst2[1]
                    # # strip gets rid of the spaces before and after
                    # pes_label_dct[spc.strip()] = label

        if 'Bimolecular' in line:

            line_lst = line.strip().split()

            if line_lst[0].strip() == 'Bimolecular' and '!' not in line:
                # Get label
                label = line_lst[1]

                # Get energy
                for line2 in input_lines[idx:]:
                    if 'Dummy' in line2:
                        ene = -10.0
                        break
                    if 'GroundEnergy' in line2:
                        ene = float(line2.split()[-1])
                        break

                # Add value to dct
                energy_dct[label] = ene

                # Add value to PES dct
                # prior_line = input_lines[idx-1]
                # line_lst2 = prior_line.split('!')
                # spc = line_lst2[1]
                # # strip gets rid of the spaces before and after
                # pes_label_dct[spc.strip()] = label

        if 'Barrier' in line:

            line_lst = line.strip().split()
            if line_lst[0].strip() == 'Barrier' and '!' not in line:
                # Get label
                [tslabel, rlabel, plabel] = line_lst[1:4]

                if ('F' not in tslabel) or ('F' in tslabel and read_fake):
                    # Get energy
                    for line2 in input_lines[idx:]:
                        if 'ZeroEnergy' in line2:
                            ene = float(line2.split()[-1])
                            break

                    # Add value to dct
                    energy_dct[tslabel] = ene

                    # Amend fake labels (may be wrong)
                    if not read_fake:
                        rlabel = rlabel.replace('F', 'P')
                        plabel = plabel.replace('F', 'P')

                    # Add the connection to lst
                    conn_lst += ((rlabel, tslabel),)
                    conn_lst += ((tslabel, plabel),)

    return energy_dct, conn_lst, pes_label_dct
