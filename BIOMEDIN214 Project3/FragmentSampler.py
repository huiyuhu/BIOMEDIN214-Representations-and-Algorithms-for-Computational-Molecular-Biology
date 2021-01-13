"""
This file contains the main fragment sampling class, which performs a Monte Carlo simulated annealing procedure to fold a protein.
"""

from pyrosetta import *
from pyrosetta.rosetta import *
from pyrosetta.rosetta.core.scoring import *

init(extra_options='-mute all  -constant_seed')

import numpy as np
import utils
from Protein import Protein
from FragmentSet import FragmentSet
import math
import os
import random


class MCMCSampler(object):
    def __init__(self, protein, fragment_set, nfrags, kmer, anneal_rate):
        """
        TO DO: initialize necessary variables
        The score function is given to you (Rosetta centroid score function)
        """

        self.scorefxn = create_score_function('score3')

        # the protein object
        self.protein = protein
        self.nfrag = nfrags
        self.kmer = kmer  # equals to 3 or 9

        # protein_structure
        self.protein_final = protein

        self.simulation_round = 1

        # temp_anneal
        self.anneal_rate = anneal_rate
        if kmer == 9:
            self.temperature = 100
            self.temperature_end = 1
        elif kmer == 3:
            self.temperature = 1
            self.temperature_end = 0.1

        # Build the nfrags of fragment candidates for each position in the protein
        self.fragment_set = fragment_set
        self.nfrag_candidates_dict = {}
        for pos in range(1, len(self.protein.sequence) - self.kmer + 1):
            self.nfrag_candidates_dict[pos] = self.fragment_set.get_lowRMS_fragments(pos, nfrags)

        # The sampled fragments indexes for each position, which is defined at self.nfrag_candidates_dict
        # self.sampled_fragments = {}
        # for pos in range(1, len(self.protein.sequence) + 1 - self.kmer):
        #    self.sampled_fragments[pos] = set()

        # For log file
        self.energy_log = self.compute_energy(self.protein)
        self.iteration_log = 0

    def compute_energy(self, protein):
        """
        TO DO
        Compute energy of protein.
        Hint: look at utils.py
        --------
        Params:
            - protein (Protein object): protein to score
        Return:
            - energy of conformation (float)
        """
        pose = protein.pose
        scorefxn = self.scorefxn
        return utils.score_pose(pose, scorefxn)

        pass

    def perturb_fragment(self, protein, frag_position):  # you may want to add more arguments
        """
        TO DO
        Sample from possible fragments for a position, and replace torsion angles of that fragment in the protein.
        ---------
        Params:
            protein: a protein project
            frag_position: the protein position that start perturbation
        Returns:
            The protein object after replacing torsion angels of the fragment
        """

        # To avoid the duplicate sampling, only sample the unsampled fragments
        # fragments_not_sampled = set(range(self.nfrag)) - self.sampled_fragments[frag_position]
        # selected_fragment = random.choice(list(fragments_not_sampled))
        # self.sampled_fragments[frag_position].add(selected_fragment)

        selected_fragment = random.choice(range(self.nfrag))

        # Conduct the perturbation on the selected fragment
        new_positions = self.nfrag_candidates_dict[frag_position][selected_fragment]
        protein_torsion = Protein(pose=protein.pose)
        for i in range(len(new_positions)):
            protein_torsion.set_torsion((frag_position + i), new_positions[i][0], new_positions[i][1])

        return protein_torsion

    def metropolis_accept(self, protein_before, protein_after, temp):  # you may want to add more arguments
        """
        TO DO
        Calculate probability of accepting or rejecting move based on Metropolis criterion.
        --------
        Params:
            protein_before: protein before the step
            protein_after: protein after the step
        Returns:
            A probability of acceptance (P_accept)
        """
        # Calculate the change in energy produced by the move
        energy_after = self.compute_energy(protein_after)
        energy_before = self.compute_energy(protein_before)
        delta_energy = energy_after - energy_before

        if delta_energy <= 0:
            prob_accept = 1
        else:
            prob_accept = math.exp(-delta_energy / temp)  # k = 1 for this project
        return prob_accept

    def anneal_temp(self, temp):
        """
        TO DO
        Anneal temperature using exponential annealing schedule. Consider kT to be a single variable (i.e. ignore Boltzmann constant)
        --------
        Params:
            temp: temperature before simulated annealing (T start)
        Returns:
            temp_anneal: temperature after simulated annealing
        """

        temp_anneal = self.anneal_rate * temp
        return temp_anneal

    def step_helper(self):
        """
        Perform the single step (single loop)
        :return:
        """
        # step 1: sample position in chain
        # sample a k-mer
        protein_seq = self.protein.sequence
        position_list = []

        # Filter out the position with nfrags samples
        # for pos in range(1, len(protein_seq) - self.kmer):  # sample position +1 ???
        #    if len(self.sampled_fragments[pos]) != self.nfrag:
        #        position_list.append(pos)
        # if len(position_list) == 0:
        #    return None, None, None
        # position_sample = random.choice(position_list)
        position_sample = random.choice(range(1, len(protein_seq) - self.kmer + 1))

        # step 2: sample fragment at position sample and replace torsions in a copied protein
        protein_torsion = self.perturb_fragment(self.protein, position_sample)

        # step 3: measure energy after replacing fragment
        energy_torsion = self.compute_energy(protein_torsion)

        # step 4: use Metropolis
        prob_accept = self.metropolis_accept(self.protein, protein_torsion, self.temperature)

        return protein_torsion, energy_torsion, prob_accept

    def step(self):
        """
        TO DO
        Take a single MCMC step. Each step should do the following:
        1. sample position in chain
            - Note: think about positions you can sample a k-mer fragment from. 
              For example, you cannot sample from position 1 because there is no phi angle
        2. sample fragment at that position and replace torsions in a *copied version* of the protein
        3. measure energy after replacing fragment
        4. accept or reject based on Metropolis criterion
            - if accept: incorporate proposed insertion and anneal temperature
            - if reject: sample new fragment (go to step 3)
        """

        iteration = 1
        while True:
            protein_torsion, energy_torsion, prob_accept = self.step_helper()
            # print("Simulation Round", self.simulation_round, iteration, energy_torsion, self.temperature)
            iteration += 1
            if protein_torsion is None:
                return
            prob_random = np.random.uniform(0, 1)
            temperature_current = self.temperature

            # keep track of the best (lowest-energy) structure produced during this simulation
            if energy_torsion < self.compute_energy(self.protein_final):
                self.protein_final = protein_torsion

            if prob_accept > prob_random:
                self.protein = protein_torsion
                self.temperature = self.anneal_temp(temperature_current)
                self.energy_log = energy_torsion
                self.iteration_log += 1
                break

    def simulate(self):
        """
        TO DO
        Run full MCMC simulation from start_temp to end_temp. 
        Be sure to save the best (lowest-energy) structure, so you can access it after.
        It is also a good idea to track certain variables during the simulation (temp, energy, and more).
        -------- 
        Params:
            none
        Returns:
            logs
        """
        temperature_all = []
        energy_all = []
        iteration_all = []
        while self.temperature >= self.temperature_end:
            self.step()
            temperature_all.append(self.temperature)
            energy_all.append(self.energy_log)
            iteration_all.append(self.iteration_log)
            self.simulation_round += 1

        return temperature_all, energy_all, iteration_all


