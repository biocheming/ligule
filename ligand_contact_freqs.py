""" Samuel D. Lotz - Lab of Dr. Alex Dickson - Michigan State
University 2015-2-16

 This script takes in a pdb, trajectories, a cutoff distance (in nm), and a
ligand name and computes the number of times in a trajectory that the
ligand and each residue are within the cutoff distance to each other.

Writes a csv of 3 columns: index of the protein residue, index of the
ligand "residue", and number of contacts.
"""

import itertools as it
import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
import sys
import pandas as pd

# the file the final data will be output to as csv format
output_f = sys.argv[1]
# the resName letter code in the PDB
ligand_code = sys.argv[2]
# the cutoff distance for accepting a ligand as in contact
cutoff = float(sys.argv[3])
# the path to the topology file
topo_path = sys.argv[4]
# the number of frames to count contacts for i.e. traj[0:frames]
# pass 'all' if to use all frames
frames = sys.argv[5]
# the trajectories to count for
trajs = [dcd for dcd in sys.argv[6:]]

# define the non-protein heteroatoms
hetero_atoms = ['HOH', 'CLA', 'CAL', ligand_code]

# load the topology
topo = md.load(topo_path).top
# load into dataframe
traj_df, traj_bonds = topo.to_dataframe()
# get the dataframe of the ligand residue
lig_df = traj_df[traj_df['resName'] == ligand_code]
# get the unique residue indices
lig_res_ix = np.unique(np.array(lig_df['resSeq']))
# get the dataframe of the protein residues
prot_df = traj_df[~traj_df['resName'].isin(hetero_atoms)]
# get the unique residue indices
prot_res_ix = np.unique(np.array(prot_df['resSeq']))
# make the product (pairwise combinations of each iterable of indices
pairs_ix = np.array([[i[0], i[1]] for i in it.product(prot_res_ix, lig_res_ix)])
# dataframe for tabulating number of times a pair is within the cutoff, a contact
pairs_df = pd.DataFrame(pairs_ix)
# add a column for number of contacts below cutoff and set to zero
pairs_df['contacts'] = 0

# go through each trajectory and count up total contacts for all of
# them
for traj in trajs:
    # load the trajectory file using the associated topology file
    traj = md.load(traj, top=topo)
    # use only the specified frames
    if frames != 'all':
        traj = traj[0:frames]
    # then compute contacts using the residue indices
    contact_dist, p = md.compute_contacts(traj, contacts=pairs_ix, scheme="closest-heavy")
    del p
    # go through these and sum up how many were below the cutoff
    num_contacts = []
    for col_i in range(contact_dist.shape[1]):
        col = contact_dist[:,col_i]
        num_contacts.append(np.where(col < cutoff)[0].shape[0])
    # add the new contacts to the pairs DataFrame
    pairs_df['contacts'] = pairs_df['contacts'] + num_contacts

# write the dataframe to a file
contacts = pairs_df['contacts']
contacts.index.name = "Residue Index"
contacts.name = ligand_name
contacts.to_csv(output_f, index=True)
