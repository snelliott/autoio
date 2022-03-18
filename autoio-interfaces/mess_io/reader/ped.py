"""
    Read the product energy distribution and store the distributions obtained
"""

import sys
from unicodedata import name
import numpy as np
import pandas as pd
import autoparse.find as apf
import autoparse.pattern as app
from ioformat import remove_comment_lines
from mess_io.reader._label import name_label_dct
from automol.util.dict_ import invert

def ped_names(input_str):
    """ Reads the ped_species and ped_output from the MESS input file string
        that were used in the master-equation calculation.
        :param input_str: string of lines of MESS input file
        :type input_str: str
        :return ped_species: list of names of connected REACS/PRODS
        :rtype: list(list(str))
        :return ped_output: output file with the product energy distribution
        :rtype: str
    """

    # Get the MESS input lines
    input_str = remove_comment_lines(input_str, delim_pattern=app.escape('!'))
    mess_lines = input_str.splitlines()
    check = 0
    ped_species = ()
    ped_output = None
    for line in mess_lines:
        if 'PEDSpecies' in line:
            check += 1
            ped_species_all = line.strip().split()[1:]
            for coupled_species in ped_species_all:
                ped_species += (tuple(coupled_species.split('_')),)
        if 'PEDOutput' in line:
            check += 1
            ped_output = line.strip().split()[1]
        if check == 2:
            # found everything you need: exit from loop
            break

    if not ped_species or not ped_output:
        print('Error: PEDSpecies and PEDOutput options incomplete. '
              'Exiting now')
        sys.exit()

    return ped_species, ped_output


def get_ped(pedoutput_str, ped_spc, energy_dct, sp_labels='inp'):
    """ Read `PEDOutput` file and extract product energy distribution at T,P.
        Energy in output set with respect to the ground energy of the products

        :param pedoutput_str: string of lines of ped_output file
        :type pedoutput_str: str
        :param pedspecies: species of interest in pedoutput
        :type pedspecies: list(list(str))
        :param energy_dct: energies of ped PES
        :type energy_dct: {label: energy} (str)
        :param sp_labels: type of pedspecies labels: 'inp' is how you find them
                in mess input, 'out' is how they are labeled in the output
        :type sp_labels: str
        :return ped_df_dct: dct(dataframe(columns:P, rows:T))
                            with the Series of energy distrib
        :rtype ped_df_dct: {((reacs,),(prods,),(None,)): dataframe(series(float))}
    """

    ped_lines = pedoutput_str.splitlines()

    # apf.where data of interest are
    pressure_i = apf.where_in('pressure', ped_lines)
    temperature_i = apf.where_in('temperature', ped_lines)

    # get T, P list
    pressure_lst = np.array([ped_lines[P].strip().split('=')[1]
                             for P in pressure_i], dtype=float)
    temperature_lst = np.array([ped_lines[T].strip().split('=')[1]
                                for T in temperature_i], dtype=float)
    # get label dictionary
    lbl_dct = name_label_dct(pedoutput_str)
    inv_lbl_dct = invert(lbl_dct)

    ped_df_dct = dict.fromkeys(energy_dct.keys())
    # find the energy for the scaling: everything refers to the products
    for spc in ped_spc:
        # allocate empty dataframe
        ped_df = pd.DataFrame(index=list(set(temperature_lst)),
                              columns=list(set(pressure_lst)), dtype=object)
        reacs, prods = spc
        label = ((reacs,), (prods,), (None,)) # rxn name
        # relabel if necessary
        if sp_labels == 'inp':
            label_messout = '->'.join([inv_lbl_dct[reacs], inv_lbl_dct[prods]])
        elif sp_labels == 'out':
            label_messout = '->'.join(spc)
        else:
            print('*Error: sp_labels must be "inp" (as in mess input) \
                or "out" (as in mess output)')
            sys.exit()

        # 0th of the energy: products energy
        ene0 = -energy_dct[prods]

        species_i = apf.where_in(label_messout, ped_lines)+1
        empty_i = apf.where_is('', ped_lines)
        final_i = np.array([empty_i[i < empty_i][0] for i in species_i])

        # column label
        column_i = apf.where_is(
            label_messout, ped_lines[species_i[0]-1].strip().split()[1:])[0]
        # reset species_i and final_i based on correspondence with label

        # check that length of pressure list is the same as new species_i
        # extract the data
        for i in np.arange(0, len(species_i)):
            # if labels don't match: go to next loop
            if label_messout not in ped_lines[species_i[i]-1]:
                continue

            pressure, temp = pressure_lst[i], temperature_lst[i]
            i_in, i_fin = species_i[i], final_i[i]

            en_prob_all = np.array(
                [line.strip().split() for line in ped_lines[i_in:i_fin]],
                dtype=float).T
            energy = en_prob_all[:][0] + ene0
            probability = en_prob_all[:][column_i]
            # build the series and put in dataframe after
            # removing negative probs and renormalizing
            # integrate with the trapezoidal rule
            prob_en = pd.Series(probability, index=energy, dtype=float)

            if len(prob_en[prob_en < 0]) > 0:
                # if there are negative values of the probability: remove them
                prob_en = prob_en[:prob_en[prob_en < 0].index[0]]

            prob_en = prob_en.sort_index()
            # integrate with trapz
            norm_factor = np.trapz(prob_en.values, x=prob_en.index)
            ped_df.loc[temp][pressure] = prob_en/norm_factor

        ped_df_dct[label] = ped_df

    return ped_df_dct
