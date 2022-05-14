"""
    Read the product energy distribution and store the distributions obtained
"""

import sys
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
                ped_species += (tuple(coupled_species.split('=')),)
        if 'PEDOutput' in line:
            check += 1
            ped_output = line.strip().split()[1]
        if check == 2:
            # found everything you need: exit from loop
            break

    if not ped_species or not ped_output:
        print('*Warning: PEDSpecies and PEDOutput options incomplete. \n'
              'Why did you call this function then? \n')

    return ped_species, ped_output


def get_ped(pedoutput_str, ped_spc, energy_dct, sp_labels='auto', hotwells=[]):
    """ Read `PEDOutput` file and extract product energy distribution at T,P.
        Energy in output set with respect to the ground energy of the products

        :param pedoutput_str: string of lines of ped_output file
        :type pedoutput_str: str
        :param ped_spc: species of interest in pedoutput
        :type ped_spc: list(list(str))
        :param energy_dct: energies of ped PES
        :type energy_dct: {label: energy} (str)
        :param sp_labels: type of pedspecies labels: 'inp' is how you find them
                in mess input, 'out' is how they are labeled in the output,
                'auto' sets inp if it finds the labels
        :type sp_labels: str
        :param hotwells: list of hot wells to find in the output
        :type hotwells: list(str)
        :return ped_df_dct: dct(dataframe(columns:P, rows:T))
                            with the Series of energy distrib
        :rtype ped_df_dct: {((reacs,),(prods,),(None,)): dataframe(series(float))}
    """
    def prob_en_single(probability, energy):
        # build the series and put in dataframe after
        # removing negative probs and renormalizing
        # integrate with the trapezoidal rule
        prob_en = pd.Series(probability, index=energy, dtype=float)
        prob_en = prob_en.sort_index()

        # if there are still negative energies (might happen for multiple ped prods)
        # remove them
        prob_en = prob_en[prob_en.index > 0]

        if len(prob_en[prob_en < 0]) > 0:
            # if there are negative values of the probability: remove them
            prob_en = prob_en[prob_en[prob_en < 0].index[-1]+1:]
        # integrate with trapz
        prob_en /= np.trapz(prob_en.values, x=prob_en.index)
        # if issues : keep only max value like dirac delta
        if any(np.isnan(prob_en.values)) or any(np.isinf(prob_en.values)):
            prob_en = np.NaN
            
        try:
            if max(prob_en) > 1:
                prob_en /= max(prob_en)
        except ValueError:
            prob_1 = apf.where_is(1, probability)
            if len(prob_1) > 1:
                prob_1 = np.array([prob_1[0]])

            prob_en = pd.Series(probability[prob_1], index=energy[prob_1])
            if prob_en.empty or energy[prob_1] < 0:
                prob_en = np.NaN
                
        except TypeError:
            pass

        return prob_en

    def indexes(label_messout, ped_lines):
        species_i = apf.where_in(label_messout, ped_lines)+1
        empty_i = apf.where_is('', ped_lines)
        final_i = np.array([empty_i[i < empty_i][0] for i in species_i])

        return species_i, final_i

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
    if lbl_dct:
        inv_lbl_dct = invert(lbl_dct)
    if sp_labels == 'auto':
        sp_labels = 'out'*(not lbl_dct) + 'inp'*(not not lbl_dct)

    ped_df_dct = {}
    # ped_df_dct = dict.fromkeys(energy_dct.keys())
    # find the energy for the scaling: everything refers to the products
    for spc in ped_spc:
        # allocate empty dataframe
        ped_df = pd.DataFrame(index=list(set(temperature_lst)),
                              columns=list(set(pressure_lst)), dtype=object)
        reacs, prods = spc
        label = ((reacs,), (prods,), (None,))  # rxn name
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

        species_i, final_i = indexes(label_messout, ped_lines)
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

            ped_df[pressure][temp] = prob_en_single(probability, energy)

        ped_df_dct[label] = ped_df

    for hotwell in hotwells:

        # relabel if necessary
        if sp_labels == 'inp':
            label_messout = 'Initial well: {}'.format(inv_lbl_dct[hotwell])
        elif sp_labels == 'out':
            label_messout = 'Initial well: {}'.format(hotwell)
        else:
            print('*Error: sp_labels must be "inp" (as in mess input) \
                or "out" (as in mess output)')
            sys.exit()
        # find products label
        species_i, final_i = indexes(label_messout, ped_lines)

        prods_list = ped_lines[species_i[0]].split('mol')[1].split()
        if sp_labels == 'inp':
            prods_outinp = pd.Series([lbl_dct[prod]
                                      for prod in prods_list], index=prods_list)
        elif sp_labels == 'out':
            prods_outinp = pd.Series(prods_list, index=prods_list)
        ene0_all = pd.Series(-np.array([energy_dct[prods_outinp[prods]]
                                        for prods in prods_list]), index=prods_list)
        # allocate dataframes and labels
        for prods in prods_outinp.values:
            label = ((hotwell,), (prods,), (None,))
            ped_df_dct[label] = pd.DataFrame(index=list(set(temperature_lst)),
                                             columns=list(set(pressure_lst)), dtype=object)

        # extract the data
        for i in np.arange(0, len(species_i)):
            i_in, i_fin = species_i[i]+1, final_i[i]
            pressure, temp = pressure_lst[pressure_i <
                                          i_in][-1], temperature_lst[temperature_i < i_in][-1]

            init_energy = float(ped_lines[i_in -
                                          2].split('=')[-1].strip()) - energy_dct[hotwell]
            en_prob_all = np.array(
                [line.strip().split() for line in ped_lines[i_in:i_fin]],
                dtype=float).T

            for pi, prods in enumerate(prods_list):

                label = ((hotwell,), (prods_outinp[prods],), (None,))
                try:
                    ped_df_dct[label][pressure][temp][init_energy]
                except TypeError:
                    ped_df_dct[label][pressure][temp] = pd.Series(dtype=object)
                except KeyError:
                    pass
                finally:
                    energy = en_prob_all[:][0] + ene0_all[prods]
                    probability = en_prob_all[:][pi+1]

                    ped_df_dct[label][pressure][temp][init_energy] = prob_en_single(
                        probability, energy)
                    ped_df_dct[label][pressure][temp] = ped_df_dct[label][pressure][temp].dropna()

    return ped_df_dct
