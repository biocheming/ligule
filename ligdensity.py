import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
import sys

output_f = sys.argv[1]
ligand_name = sys.argv[2]
cutoff = float(sys.argv[3])
topo_path = sys.argv[4]
trajs = [dcd for dcd in sys.argv[5:]]

pairs_hist = np.array([])
bincounts = np.array([])
output = np.array([])
close_counts = []
for traj in trajs:
    # load the trajectory file using the associated topology file
    t0 = md.load(traj, top=topo_path)

    # the index of the protein atoms
    t0_prot_ix = t0.top.select_atom_indices(selection='minimal')
    # the indices of the DMSO atoms
    t0_ligand_ix = t0.top.select("resname {0}".format(ligand_name))

    # get a (2,n_pairs) array for input into compute distances
    pairs_ix = t0.top.select_pairs(selection1=t0_prot_ix, selection2=t0_ligand_ix)
    
    # calculate the distances returns (n_frames, n_pairs) array
    pair_distances = md.compute_distances(t0, pairs_ix)
    
    # get the indices of those under the cutoff
    cutoff_indices = np.where(pair_distances < cutoff)
    # go through these and sum up how many were below the cutoff
    for col in range(cutoff_indices[1].shape):
        num_close = np.where(col < cutoff)
    # returns a tuple for each axis index (frames, pair)
    pairs_hist = np.append(pairs_hist, cutoff_indices[1])
    # pair this with the index of the atom pairs
    output = np.concatenate((pairs_ix, pairs_hist), axis=1)
    #bincounts = np.append(bincounts, np.bincount(cutoff_indices[1]))

# normalize based on the longest distance
#norm_bincount = bincount/bincount.max()
# plot using our own histogram method
#left = [left for left in range(len(norm_bincount))]
#plt.bar(left, norm_bincount)

np.savetxt(output_f+".csv", output, delimiter=",")

# plot a histogram
plt.hist(pairs_hist, bins=len(t0_prot_ix), range=(0, len(t0_prot_ix)), normed=True)
plt.savefig(output_f+".png", bbox_inches='tight')
