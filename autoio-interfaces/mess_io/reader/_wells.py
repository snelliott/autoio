""" Read the merged wells from the auxiliary file
"""

import numpy
from phydat import phycon
import autoparse.pattern as app
import autoparse.find as apf


def merged_wells(mess_aux_str, pressure, temp):
    """ Parse the auxiliary MESS output file string for all of the groups
        of wells which merged at a given pressure and temperature
        from a Master Equation simulation

        :param mess_aux_str: string for the auxiliary file
        :type mess_aux_str: str
        :rtype: tuple(tuple(str))
    """

    # Convert input pressure to bar for reading
    pressure *= phycon.ATM2BAR

    # Get the lines
    mess_lines = mess_aux_str.splitlines()

    # Parses number of wells and lines they appear on for second loop
    # Done for each set of wells at each T and P
    merged_well_lines = []
    for i, line in enumerate(mess_lines):
        if 'number of species =' in line:
            cond_line = mess_lines[i-1].strip().split()
            well_pressure = float(cond_line[2])
            well_temp = float(cond_line[-2])
            if (numpy.isclose(pressure, well_pressure) and
                    numpy.isclose(temp, well_temp)):
                nspc = int(line.strip().split()[-1])
                merged_well_lines.append((nspc, i))

    # Find all the lines in the block to grab the wells
    merged_well_lst = ()
    for nspc, line_idx in merged_well_lines:
        for idx in range(nspc):
            well_line = mess_lines[line_idx+1+idx*2]
            well_names = well_line.strip().split()
            if len(well_names) > 1:
                merged_well_lst += (tuple(well_names),)

    return merged_well_lst


def well_average_energy(log_str, well, temp):
    """ Obtain the average energy of each well from the output
        of MESS rate calculations.

        Returns the energies in hartrees.

        :param log_str: string of the MESS .log file
        :type log_str: str
        :param well: name of the well to get energy for
        :type well: str
        :param temp: temperature to get the energy for
        :type temp: float
        :rtype: dict[str: float]
    """

    # Loop through file to find the energy block with requested temp
    block_ptt = (
        'MasterEquation::set:  starts' +
        app.capturing(app.one_or_more(app.WILDCARD, greedy=False)) +
        'MasterEquation::set:  done'
    )
    blocks = apf.all_captures(block_ptt, log_str)

    if blocks is not None:
        for i, block in enumerate(blocks):
            templine = block.splitlines()[1]
            blocktemp = float(templine.strip().split()[2])
            if numpy.isclose(blocktemp, temp, atol=0.01):
                block_idx = i
                break
    else:
        block_idx = None

    # If block with requested temp found, get energies
    ene_dct = None
    if block_idx is not None:
        ptt = (
            app.capturing(app.VARIABLE_NAME) + app.SPACE +
            'Well:' + app.SPACE +
            'average energy =' + app.SPACE +
            app.capturing(app.NUMBER) + app.SPACE +
            'kcal/mol'
        )
        caps = apf.all_captures(ptt, blocks[block_idx])
        ene_dct = dict(
            ((cap[0], float(cap[1])*phycon.KCAL2EH) for cap in caps)
        )
        ene = ene_dct[well]

    return ene
