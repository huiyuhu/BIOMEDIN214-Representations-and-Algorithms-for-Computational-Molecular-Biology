"""
This is the master file, which you should use to set up and run the simulations.
You may define functions or classes as necessary

For an input sequence, do the following (see project page for details):
	1. load sequence from fasta file and initialize protein into extended configuration
	2. Run a series of simulations (n=5 or 10):
		- Perform MCMC sampling with 9-mer fragments from kT=100 to kT=1 (assembly stage)
		- Perform MCMC sampling with 3-mer fragments from kT=1 to kT=0.1 (refinement stage)
		- Take best (lowest-energy) structure after refinement and perform energy minimization (see utils.relax)
		- Log energy and RMSD to native structure after minimization
	3. Visualize lowest-RMSD structure in PyMol

"""

from pyrosetta import *
from pyrosetta.rosetta import *

init(extra_options='-mute all -constant_seed')

from Protein import Protein
from FragmentSampler import MCMCSampler
from FragmentSet import FragmentSet

import utils
import argparse
import pandas as pd
import numpy as np
import time
import os


def main():
    # Read in arguments
    parser = argparse.ArgumentParser(description='Performs ab initio protein structure prediction')
    parser.add_argument('--fasta', help='.fasta file containing the sequence')
    parser.add_argument('--logdir', help='directory to save all log files', default='.')
    parser.add_argument('--nsims', type=int, help='number of simulations', default=1)
    parser.add_argument('--nfrags', type=int, help='number of fragments to sample at each iteration', default=3)
    parser.add_argument('--anneal_rate', type=float, help='temperature annealing parameter', default=0.999)
    args = parser.parse_args()
    # according to input, get the file, sequence
    fasta_file = args.fasta
    #print("FASTA!!!!!!!" + args.fasta)
    fasta = open(fasta_file)
    protein_seq = fasta.readlines()[1].strip()
    protein_name = args.fasta.split(sep='.')[0]
    # protein object
    protein = Protein(sequence=protein_seq)

    # frag & rmsd file name
    frag9mer_file = protein_name + '_9mers.frag'
    frag3mer_file = protein_name + '_3mers.frag'
    rmsd9mer_file = protein_name + '_9mers.rmsd'
    rmsd3mer_file = protein_name + '_3mers.rmsd'
    fragment_set_9mer = FragmentSet(frag9mer_file, rmsd9mer_file)
    fragment_set_3mer = FragmentSet(frag3mer_file, rmsd3mer_file)

    sim_number = []
    energy = []
    rmsd = []
    for i in range(args.nsims):
        # stage 1: 9-mer
        sampler_stg1 = MCMCSampler(protein=protein, fragment_set=fragment_set_9mer, kmer=9, nfrags=args.nfrags,
                                   anneal_rate=args.anneal_rate)
        temp_log_stg1, energy_log_stg1, iter_log_stg1 = sampler_stg1.simulate()

        protein_final_9mer = sampler_stg1.protein_final

        # stage 2: 3-mer, starting with the best structure from stage 1
        sampler_stg2 = MCMCSampler(protein=protein_final_9mer, fragment_set=fragment_set_3mer, kmer=3,
                                   nfrags=args.nfrags, anneal_rate=args.anneal_rate)
        temp_log_stg2, energy_log_stg2, iter_log_stg2 = sampler_stg2.simulate()
        protein_final_3mer = sampler_stg2.protein_final

        # save the final structure into pdb file
        protein_final_3mer.save_pdb(filename='best.pdb')

        # perform energy minimization on the output structure (relaxation)
        # - pdb (str): path to input structure in PDB format
        # - native (str): path to native structure in PDB format

        protein_final, rmsd_final, score_final = utils.relax(pdb='best.pdb', native=protein_name + '.pdb')

        protein_final.save_pdb(filename='best_fast_relax.pdb')

        # Save the results and log
        sim_number.append(i)
        energy.append(score_final)
        rmsd.append(rmsd_final)
    # print(sim_number, energy, rmsd)
    sim_log = pd.DataFrame(list(zip(sim_number, energy, rmsd)))
    sim_log.columns = ['sim_number', 'energy', 'rmsd']
    sim_log.to_csv((args.logdir + '/simulation_summary.txt'), sep='\t', index=False)


if __name__ == '__main__':
    main()
    pass
