import pandas as pd

import utils
import numpy as np
import os


class FragmentSet(object):
    def __init__(self, fragfile, rmsdfile):
        """
		This class contains the fragment library for the input protein. It must do the following:
		- Read in fragment file and parse fragments at each position. Fragment files are of the form <protein>_<frag_len>mers.frag
		- Read in RMSD file containing pre-calculated RMSD to native structure for each fragment at each position.
		- Based on fragments and their corresponding RMSDs, rank fragments at each position by RMSD
		
		"""
        ## The only important data are the phi and psi columns
        ## 19-27  -- phi; 28-36  -- psi

        rmsd_data = pd.read_csv(rmsdfile, sep='\t', header=None)
        self.rmsd = rmsd_data.sort_values(2, ascending=True).groupby(0)

        pos_dict = {}
        with open(fragfile) as my_file:
            current_pos = ""
            for line in my_file:
                line_items = line.strip().split()
                if len(line_items) == 0:
                    continue
                elif line_items[0] == 'position:':
                    current_pos = line_items[1]
                    pos_dict[current_pos] = {}
                else:
                    frag_id = line_items[-1]
                    if frag_id[0] == 'F':
                        frag_id = frag_id[1:]
                    if frag_id not in pos_dict[current_pos]:
                        pos_dict[current_pos][frag_id] = []
                    phi = float(line_items[5])
                    psi = float(line_items[6])
                    pos_dict[current_pos][frag_id].append((phi, psi))

        self.position_frag_dict = pos_dict
        pass

    def get_lowRMS_fragments(self, pos, N):
        """
		Returns the top-ranked fragments by RMSD at a defined position in the chain
		--------
		Params
			- pos (int): fragment position in chain (1-indexed)
			- N (int): number of fragments to return
		Returns
			- lowRMS_fragments (list): top N fragments at pos by RMSD. This should be a list of lists of (phi, psi) tuples. 
			  For example, a 3-mer fragment could be represented as the following: [(-60.892, 142.456), (-72.281, 128.933), (-132.337, -175.477)]
		"""
        frag_id_list = self.rmsd.head(N).loc[self.rmsd.head(N)[0] == pos][1].values.tolist()

        lowRMS_fragments = []
        for frag_id in frag_id_list:
            phi_psi = self.position_frag_dict[str(pos)][str(frag_id)]
            lowRMS_fragments.append(phi_psi)

        return lowRMS_fragments

        pass
