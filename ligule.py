from __future__ import print_function
import click

# While functions are significantly changed inspiration for functions
# from pdbtools

# reference data section
BORING_HETEROATOMS = ['HOH','CLA','CAL', 'CA', 'SO4', 'IOD', 'NA', 'CL', 'GOL', 'PO4'] # heteroatoms in the HET label which are not interesting, e.g. solvent
HETEROATOMS = ['HOH','CLA','CAL'] # heteroatoms that are labelled as atoms
# the command subgroup
@click.group()
def cli():
    pass

def heaviest_atom(atoms_df):
    """ Given a dataframe of atoms return the one with the largest mass."""
    import mdtraj as md
    heaviest_atom_ix = 0
    heaviest_atom_mass = 0.0
    for row in atoms_df.iterrows():
        element = row[1]['element']
        mass = md.core.element.Element.getBySymbol(element).mass
        if mass > heaviest_atom_mass:
            heaviest_atom_ix = row[0]
            heaviest_atom_mass = mass
    return heaviest_atom_ix


# the contact frequency command
#
@click.command()
@click.option('--output', default=None, help="file path to save contact frequency output")
@click.option('--cutoff', default=0.6, type=float, help="pairwise protein residue-ligand distances under this cutoff (in nm) counted as in contact")
@click.option('--nframes', default='all', type=str, help="list of the number of frames to take from the respective trajectories use 'all' to get all frames if only one is given with multiple trajectories it is used for all trajectories, e.g. 1500 or all or 1000,all,all, default is 'all' ")
@click.option('--hetatoms', default=','.join(HETEROATOMS), help="list of non-protein heteroatom codes not including the ligand labelled as ATOM in the pdb that should be HET, e.g. 'HOH,CLA', default={0}".format(','.join(HETEROATOMS)))
@click.option('--plot', is_flag=True, help="produce a matplotlib bar chart of the frequencies")
@click.option('--weights', default=None, type=str, help="paths to the weights files e.g. /path/1.dat,/path/2.dat")
@click.option('--byRes', 'by_particle', flag_value="residue", default=True, help="output contacts on a residue basis")
@click.option('--byAtom', 'by_particle', flag_value="atom", default=False, help="output contacts on an atom basis")
@click.option('--ligAtom', default=None, type=int, help="if byAtom was enabled choose the atom by PDB index to get distances from the protein")
@click.argument('ligand', nargs=1, type=str)
@click.argument('topology', nargs=1, type=click.Path(exists=True))
@click.argument('traj_paths', nargs=-1, type=click.Path(exists=True))
def contact_freqs(output, cutoff, nframes, hetatoms, plot, weights, by_particle, ligatom, ligand, topology, traj_paths):
    import mdtraj as md
    import numpy as np
    import pandas as pd
    import itertools as it
    import sys
    # rename parameters
    ligand_code = str(ligand)
    del ligand
    topo_path = topology
    del topology
    weight_paths = weights
    del weights
    
    # make sure the byAtom flag was given if the ligAtom was given
    try:
        if not ligatom is None and not by_particle == 'atom':
            raise ValueError
    except ValueError:
        click.echo("Ligand atom given without byAtom flag, ignoring")

    # parse the weights option if given
    if weight_paths:
        tmp = [weight_path for weight_path in weight_paths.split(',')]
        weight_paths = tmp
        del tmp
    # if the 
    # load the number of frames list and convert integer string to
    # ints unless its value is 'all'
    tmp = []
    for nframe in nframes.split(','):
        try:
            tmp.append(int(nframe))
        except ValueError:
            if nframe == 'all':
                tmp.append(nframe)
            else:
                raise
    nframes = tmp
    del tmp

    # expand frames slice to others if only one is given
    if len(nframes) == 1:
        nframes = [nframes[0] for f in range(len(traj_paths))]
    # check to make sure that the correct number of frames were given
    if len(nframes) != len(traj_paths):
        raise ValueError("Wrong number of frame slices for trajectories given")

    # read in the heteroatoms to a list
    hetero_atoms = [str(het) for het in hetatoms.split(',')]
    # add the ligand to the heteroatoms
    hetero_atoms.append(ligand_code)
    
    # load the topology
    topo = md.load(topo_path).top
    # load into dataframe
    traj_df, traj_bonds = topo.to_dataframe()
    # get the dataframe of the ligand residue if it is in the topology
    try:
        if ligand_code not in traj_df['resName'].values:
            raise ValueError
    except ValueError:
        sys.exit("supplied ligand code, '{0}', is not in the topology".format(ligand_code))

    # make pairs of particle types from the protein based on the
    # 'by_particle' flag
    if by_particle == 'residue':
        # column named 'resSeq', matches each residue from protein
        by_df_col = 'resSeq'
    elif by_particle == 'atom':
        # column named 'serial' for the original pdb numbering,
        # matches each protein atom
        by_df_col = 'serial'

    # get the rows that correspond to the ligand code given
    lig_df = traj_df[traj_df['resName'] == ligand_code]

    # get the dataframe of the protein particles
    prot_df = traj_df[~traj_df['resName'].isin(hetero_atoms)]

    # get the index of particles for the ligand
    if by_particle == 'residue':
        lig_ix = np.unique(np.array(lig_df['resSeq']))        
        # get the index of residues for the protein
        prot_ix = np.unique(np.array(prot_df['resSeq']))
    # if by atom is chosen
    elif by_particle == 'atom':
        # if chosen atom is not given use the heaviest atom
        if not ligatom is None:
            # compute distances to
            lig_ix = np.array([heaviest_atom(lig_df)])
        # use the given atom index
        else:
            # minus 1, mdtraj is 0 indexed and pdb is 1 indexed
            lig_ix = ligatom-1
            click.echo(topo.atom(lig_ix))
        # get the index of atoms for the protein
        prot_ix = np.array(prot_df.index)

    # make the product (pairwise combinations of each iterable of indices)
    pairs_ix = np.array([[i[0], i[1]] for i in it.product(prot_ix, lig_ix)])
    # dataframe for tabulating number of times a pair is within the cutoff, a contact
    pairs_df = pd.DataFrame(pairs_ix)
    # add a column for number of contacts below cutoff and set to zero
    pairs_df['contacts'] = 0        

        
    # go through each trajectory and count the contacts 
    for i, traj_name in enumerate(traj_paths):
        # load the trajectory file using the associated topology file
        traj = md.load(traj_name, top=topo)
        # use only the specified frames
        if nframes[i] != 'all':
            try:
                traj = traj[0:nframes[i]]
            except IndexError:
                click.echo("Trajectory has fewer frames than queried, using all frames")

        if by_particle == 'residue':
            # then compute contacts using the residue indices
            contact_dist, p = md.compute_contacts(traj, contacts=pairs_ix, scheme="closest-heavy")
            del p
        elif by_particle == 'atom':
            # then compute contacts from the atoms
            contact_dist = md.compute_distances(traj, pairs_ix)
            
        # go through these and sum up how many were below the cutoff
        num_contacts = []
        for col_i in range(contact_dist.shape[1]):
            col = contact_dist[:,col_i]
            num_contacts.append(np.where(col < cutoff)[0].shape[0])
        # add the new contacts to the pairs DataFrame values
        pairs_df['contacts'] = pairs_df['contacts'] + num_contacts

        # apply weights if supplied
        if weight_paths:
            # load weights file from option
            weights = pd.read_csv(weight_paths[i], header=None)
            # multiply counts of contacts by weights
            pairs_df['contacts'] = pairs_df['contacts'] * weights[0]

    # set contacts dataframe metadata
    contacts = pairs_df['contacts']
    contacts.index.name = "Residue Index"
    contacts.name = ligand_code

    # try to make output if desired
    if output:
        try:
            contacts.to_csv(output, index=True, index_label="Residue Index", header=True)
        except IOError as e:
            click.echo("Bad output path {0}".format(e))

    # if not send to stdout
    if output is None:
        click.echo("Not writing output")
        sys.stdout.write(contacts.to_string())
        sys.stdout.write("\n")
    
    # plotting
    if plot == True:
        click.echo("Plotting not implemented yet")

    return contacts




def pdbHetatomCodes(pdb_lines):
    """ Returns all codes for HETATOMS in a pdb file
    """
    hetatoms = list(set([l[7:10].strip() for l in pdb_lines if l.startswith("HET   ")]))

    return hetatoms

# command to list all ligands in a pdb file
@click.command()
@click.option('--ignore', default=None, help="specify a heteroatom to skip")
@click.option('--boring/--no-boring', default=False, help="include boring heteroatoms that are probably not the ligand. {0!r}".format(BORING_HETEROATOMS))
@click.argument('files', nargs=-1, type=click.Path(exists=True))
def ligs(ignore, boring, files):
    """
    Function to call if run from command line.
    """
    outs = []
    for i, pdb_file in enumerate(files):
        with open(pdb_file,'r') as f:
            pdb = f.readlines()

        pdb_id = pdb_file.split(".pdb")[0]

        hetatoms = pdbHetatomCodes(pdb)
        if not boring:
            ligands = [het for het in hetatoms if het not in BORING_HETEROATOMS]
        if boring:
            ligands = hetatoms

        ligand_out = "\n".join(ligands)
        out = "pdb {0}:\n{1}".format(pdb_id, ligand_out)
        outs.append(out)
        click.echo(out)

    return outs


# set groupings
cli.add_command(contact_freqs)
cli.add_command(ligs)

if __name__ == "__main__":
    cli()
