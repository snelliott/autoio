"""
    Read the HotEnergies distribution and store it
"""

import numpy as np
import pandas as pd
import autoparse.find as apf
from mess_io.reader._pes import pes

def get_hot_species(input_str):
    """ Reads the HotSpecies from the MESS input file string
        that were used in the master-equation calculation.
        returns dictionary with associated energies
        :param input_str: string of lines of MESS input file
        :type input_str: str
        :return hotspecies: dictionary of hotspecies energies
        :rtype: dct{hotspecies: en}
    """

    # Get the MESS input lines
    energy_dct, _, _, _ = pes(input_str)
    mess_lines = input_str.splitlines()
    hotsp_i = apf.where_in('HotEnergies', mess_lines)[0]
    num_hotsp = int(mess_lines[hotsp_i].strip().split()[1])
    #hotspecies = [None]*num_hotsp
    hotspecies_en = {}
    for line in mess_lines[hotsp_i+1:hotsp_i+1+num_hotsp]:
        hotname= line.strip().split()[0]
        hotspecies_en[hotname] = energy_dct[hotname]

    return hotspecies_en


def extract_hot_branching(hotenergies_str, hotspecies_en, species_lst,
                          temps, pressures):
    """ Extract hot branching fractions for a single species
        :param hotenergies_str: string of mess log file
        :type hotenergies_str: str
        :param hotspecies_en: dct of hotspecies and corresponding energy
        :type hotspecies_en: dct{hotspecies: en}
        :param species_lst: list of all species on the PES
        :type species_lst: list
        :return hoten_dct: hot branching fractions for hotspecies
        :rtype hoten_dct: dct{hotspecies: df[P][T]:df[allspecies][energies]}
    """
    lines = hotenergies_str.splitlines()
    # for each species: dataframe of dataframes BF[Ti][pi]
    # each of them has BF[energy][species]
    # preallocations
    hotspecies_lst = list(hotspecies_en.keys())
    hoten_dct = {s: pd.DataFrame(index=temps, columns=pressures)
                 for s in hotspecies_lst}

    # 1. for each P,T: extract the block
    pt_i_array = apf.where_in(['Pressure', 'Temperature'], lines)
    hot_i_array = apf.where_in(['Hot distribution branching ratios'], lines)
    end_hot_i_array = apf.where_in(
        ['prompt', 'isomerization', 'dissociation'], lines)

    # 2. find Hot distribution branching ratios:
    for i, hot_i in enumerate(hot_i_array):

        # extract block, PT, and species for which BF is assigned
        lines_block = lines[hot_i+2:end_hot_i_array[i]]

        _press, _temp = [
            float(var)
            for var in lines[pt_i_array[i]].strip().split()[2:7:4]]

        species_bf_i = lines[hot_i+1].strip().split()[3:]

        # for each hotspecies: read BFs
        # rescale energy by the hotspecies energy on the PES!!
        # ref 0 energy is the ref for the PES, even for the hotspecies
        for hotspecies in hotspecies_lst:
            hot_e_lvl, branch_ratio = [], []
            sp_i = apf.where_in(hotspecies, species_bf_i)
            ref_en = hotspecies_en[hotspecies]

            for line in lines_block:
                line = line.strip()
                if line.startswith(hotspecies):
                    hot_e = float(line.split()[1]) - ref_en
                    # pick only values above 0
                    if hot_e not in hot_e_lvl and hot_e > 0:
                        branch_ratio_arr = np.array(
                            list(line.split()[2:]), dtype=float)

                        # check that value of reactant branching is between 0 and 1
                        # if > 1, keep it so you can account for that anyway
                        if sp_i.size > 0:
                            if branch_ratio_arr[sp_i] < 0 or branch_ratio_arr[sp_i] > 1:
                                continue

                        # remove negative values or values >1
                        _arr = [abs(x*int(1e-5 < x <= 1))
                                for x in branch_ratio_arr]
                        br_filter = np.array(_arr, dtype=float)

                        # if all invalid: do not save
                        if all(br_filter == 0):
                            continue
                        br_renorm = br_filter/np.sum(br_filter)
                        # append values
                        branch_ratio.append(br_renorm)
                        hot_e_lvl.append(hot_e)

            hot_e_lvl = np.array(hot_e_lvl)
            branch_ratio = np.array(branch_ratio)

            # 3. allocate in the dataframe
            bf_hotspecies = pd.DataFrame(
                0, index=hot_e_lvl, columns=species_lst)
            bf_hotspecies[species_bf_i] = branch_ratio
            hoten_dct[hotspecies][_press][_temp] = bf_hotspecies

    return hoten_dct
