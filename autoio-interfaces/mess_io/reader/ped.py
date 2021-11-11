"""
    Read the product energy distribution and store the distributions obtained
"""
import sys
import numpy as np
import pandas as pd
from autoparse import find
import autoparse.pattern as app
from ioformat import remove_comment_lines

def ped_names(input_str):
    """ Reads the PEDSpecies and PEDOutput from the MESS input file string
        that were used in the master-equation calculation.
        :param input_str: string of lines of MESS input file
        :type input_str: str
        :return PEDSpecies: list of names of connected REACS/PRODS
        :rtype: list(list(str))
        :return PEDOutput: output file with the product energy distribution
        :rtype: str
    """

    # Get the MESS input lines
    input_str = remove_comment_lines(input_str, delim_pattern=app.escape('!'))
    mess_lines = input_str.splitlines()
    check = 0
    PEDSpecies = []
    PEDOutput = None
    for line in mess_lines:
        if 'PEDSpecies' in line:
            check += 1
            PEDSpecies_all = line.strip().split()[1:]
            for coupled_species in PEDSpecies_all:
                PEDSpecies.append(coupled_species.split('_'))
        if 'PEDOutput' in line:
            check += 1
            PEDOutput = line.strip().split()[1]
        if check == 2:
            # found everything you need: exit from loop
            break

    if not PEDSpecies or not PEDOutput:
        print('Error: PEDSpecies and PEDOutput options incomplete. Exiting now')
        sys.exit()

    return PEDSpecies, PEDOutput


def get_ped(pedoutput_str, pedspecies, energy_dct):
    """ Read the pedoutput file and extracts product energy distribution at T,P
        Energy in the output is set with respect to the ground energy of the products
        :param pedoutput_str: string of lines of PEDOutput file
        :type pedoutput_str: str
        :param species: species of interest in pedoutput
        :type species: list(list(str))
        :param energies: dictionary of energies of the species/barriers in the PES
        :type energies: dct(label: float)
        :param barriers_dct: barriers associated to pedspecies
        :type barriers_dct: {'REAC->PROD': barrierame} (str)

        :return ped_df_dct: dct(dataframe(columns:P, rows:T)) with the Series of energy distrib
        :rtype ped_df_dct: {'REAC->PROD': dataframe(series(float))}
    """

    ped_lines = pedoutput_str.splitlines()

    # find where data of interest are
    pressure_i = find.where_in('pressure', ped_lines)
    temperature_i = find.where_in('temperature', ped_lines)

    # get T, P list
    pressure_lst = np.array([ped_lines[P].strip().split('=')[1]
                             for P in pressure_i], dtype=float)
    temperature_lst = np.array([ped_lines[T].strip().split('=')[1]
                                for T in temperature_i], dtype=float)

    ped_df_dct = dict.fromkeys(energy_dct.keys())
    # find the energy for the scaling: everything refers to the products
    for species in pedspecies:
        # allocate empty dataframe
        ped_df = pd.DataFrame(index=list(set(temperature_lst)),
                              columns=list(set(pressure_lst)), dtype=object)
        _, prods = species
        label = '->'.join(species)
        # 0th of the energy: products energy
        E0 = -energy_dct[prods]

        species_i = find.where_in(label, ped_lines)+1
        empty_i = find.where_is('', ped_lines)
        final_i = np.array([empty_i[i < empty_i][0] for i in species_i])

        # column label
        column_i = find.where_is(
            label, ped_lines[species_i[0]-1].strip().split()[1:])[0]
        # reset species_i and final_i based on correspondence with label

        # check that length of pressure list is the same as new species_i
        # extract the data
        for i in np.arange(0, len(species_i)):
            # if labels don't match: go to next loop
            if label not in ped_lines[species_i[i]-1]:
               continue

            P, T, i_in, i_fin = [pressure_lst[i],
                                 temperature_lst[i], species_i[i], final_i[i]]

            en_prob_all = np.array(
                [line.strip().split() for line in ped_lines[i_in:i_fin]], dtype=float).T
            energy = en_prob_all[:][0] + E0
            probability = en_prob_all[:][column_i]
            # build the series and put in dataframe after removing negative probs and renormalizing
            # integrate with the trapezoidal rule
            prob_en = pd.Series(probability, index=energy, dtype=float)

            if len(prob_en[prob_en < 0]) > 0:
                # if there are negative values of the probability: remove them
                prob_en = prob_en[:prob_en[prob_en < 0].index[0]]

            prob_en = prob_en.sort_index()
            # integrate with trapz
            norm_factor = np.trapz(prob_en.values, x=prob_en.index)
            ped_df.loc[T][P] = prob_en/norm_factor

        ped_df_dct[label] = ped_df

    return ped_df_dct
