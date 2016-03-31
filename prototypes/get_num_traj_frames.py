from __future__ import print_function

""" Samuel D. Lotz - Lab of Dr. Alex Dickson - Michigan State
University 2015-2-16



"""

import mdtraj as md
import sys

# the file the final data will be output to as csv format
topo_path = sys.argv[1]

# the trajectories to count for
traj_names = sys.argv[2:]

# load the topology
topo = md.load(topo_path).top

# print the number of frames for a list of trajectories
for traj_name in traj_names:
    print(traj_name, "frames:", md.load(traj_name, top=topo).n_frames)
    
