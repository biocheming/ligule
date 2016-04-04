from __future__ import print_function

import sys
sys.path.insert(0, "/usr/lib/python2.7/dist-packages")


# standard library imports
import os
import os.path as osp

# internal imports
#import plip.ligcomplex as ligclx
#import plip.report as ligrep

# external imports
import click
import pandas as pd
import biopandas.pdb as pdPDB

# While functions are significantly changed inspiration for functions
# from pdbtools

# reference data section
BORING_HETEROATOMS = ['HOH','CLA','CAL', 'CA', 'SO4', 'IOD', 'NA', 'CL', 'GOL', 'PO4'] # heteroatoms in the HET label which are not interesting, e.g. solvent
HETEROATOMS = ['HOH','CLA','CAL'] # heteroatoms that are labelled as atoms

# types

FREQ_SUMS = 'freq_sums'
WEIGHTED_SUMS = 'weighted_sums'

FREQ_OUTPUT_STR_DICT = {FREQ_SUMS : 'freqSum',
                   WEIGHTED_SUMS : 'freqWeight'
                   }


# the command subgroup
@click.group()
def cli():
    pass

@click.group()
def vis():
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
    List all ligands in a pdb file, with presets.
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


@click.command()
@click.option('--clust-format/--no-clust-format', default=True, help="parses the filename for the pdb id 'lig_clust###.pdb'")
@click.option('--interaction', '-i', multiple=True, default=None, type=str, help="set which interactions you want to profile: 'hbonds': hydrogen bonds, 'hydros': hydrophobic contacts.")
@click.option('--output-dir', '-d', default=None, type=str, help="an output directory to put output files, will make new directory if it doesn't exist, but won't overwrite")
@click.option('--output-base', '-O', default=None, type=str, help="base filename to ouput results to")
@click.option('--output-type', '-T', multiple=True, default=None, type=str, help="type of output requested: 'lig_int_freq': interaction type frequency by ligand atom")
@click.argument('ligand', nargs=1, type=str)
@click.argument('files', nargs=-1, type=click.Path(exists=True))
def lig_interactions(clust_format, interaction, output_dir, output_base, output_type, ligand, files):
    """Outputs summarized data for the type of interactions a ligand has with the
    protein for many pdbs, depends on PLIP module.
    """

    if interaction is None:
        raise ValueError("You must enter at least one interaction type")

    if interaction and not output_base:
        raise ValueError("You must give an output base name.")

    # make an output dir if requested
    if output_dir:
        os.mkdir(output_dir)
        

    # DEBUG
    print("Starting interaction sets")

    outputs = None
    for pdb_path in files:

        # DEBUG
        print("Int set starting:", pdb_path)

        # Hack to set the id for the cluster
        if clust_format:
            complex_id = int(osp.splitext(osp.basename(pdb_path))[0].split("clust")[-1])
        else:
            raise ValueError("Only the cluster format 'lig_clust#.pdb' is supported currently please enter files with this name.")
            
        # create a ligand-protein complex
        lig_comp = ligclx.create_lig_complex(pdb_path, ligand)

        # get the interaction set we are interested in
        int_set = lig_comp.interaction_sets[lig_comp.interaction_sets.keys()[0]]

        # generate data row
        new_out = ligrep.complex_outputs(int_set, complex_id, interaction, output_type)

        # if it is not the first row add this to the output
        if outputs is None:
            outputs = new_out
        else:
            outputs = ligrep.concat_complex_outputs(outputs, new_out)

        # DEBUG
        print("Int set ending:", pdb_path)

    # write out the output files
    ligrep.write_complexes_output(outputs, output_dir, output_base)

    return outputs


@click.command()
@click.option('--fformat', '-f', default='pdb', help="format of ligand molecule input files")
@click.option('--output-dir', '-d', default=None, type=str, help="an output directory to put output files, will make new directory if it doesn't exist, but won't overwrite")
@click.option('--output-base', '-O', default=None, type=str, help="base filename to ouput results to")
@click.option('--vis', '-Z', default=None, multiple=True, type=(str, str), help="produce a visualization of type VIS-TYPE with the data DATA-PATH; '-Z VIS-TYPE DATA-PATH")
@click.argument('ligand-file', nargs=1, type=click.Path(exists=True))
@click.argument('freq-data', nargs=1, type=click.Path(exists=True))
def freq(fformat, output_dir, output_base, vis, ligand_file, freq_data):
    """Output (pdb) files with values set for visualization of
interaction frequency.

    Available visualization types (vis):
        '{0}' : None -- sums the frequencies
        '{1}' : [weights] -- weight frequency by the weights and sum

    """.format(FREQ_SUMS, WEIGHTED_SUMS)

    if output_base is None:
        raise ValueError("Must provide a basename")
    if vis is None:
        raise ValueError("Must provide a vis type")
    # make an output dir if requested
    if output_dir:
        os.mkdir(output_dir)

    # load the frequency data as a DataFrame
    freq_data = pd.read_csv(freq_data, index_col=0)
    print("FREQ_DATA\n", freq_data.shape)
    for vis_type, data in vis:

        # make requested visualization datas
        if FREQ_SUMS == vis_type:
            # no data to transform freq_data
            vis_data = freq_data.sum()
            
        elif WEIGHTED_SUMS == vis_type:
            # multiply each row by that complexes weight
            # TODO HACK because my input files are not standard
            weights = pd.Series.from_csv(data, index_col=1, sep=' ')
            print("WEIGHTS\n", weights.shape)
            vis_data = freq_data.rmul(weights, axis='index')
            print("VIS_DATA_WITH_WEIGHTS\n", vis_data)            
            # sum by columns for a weighted sum
            vis_data = vis_data.sum()
            print("AFTER SUM\n", vis_data)

        # output results
        if fformat == 'pdb':
            # write file with vis_type string
            outpath = "{0}_{1}.pdb".format(output_base,
                                                    FREQ_OUTPUT_STR_DICT[vis_type])

            if output_dir:                
                outpath = osp.join(output_dir, filename)

            # TODO biopandas might not update the pdb_text after you
            # alter the df so I will output pdb from the pdb_replace_bfactor file for now

            # get the text with the b_factor column
            # replaced
            pdb_text = pdb_replace_bfactor(ligand_file, vis_data, outpath)


            # with open(outpath, 'w') as fw:
            #     fw.write(pdb_text)

            

    return None

def pdb_replace_bfactor(pdb_path, series, outpath):
    """Use biopandas to read in as a dataframe and return the pdb text."""

    
    pl = pdPDB.PandasPDB().read_pdb(pdb_path)
    print("REPLACING WITH\n", series)
    print("b_factor\n",  pl.df['ATOM']['b_factor'])
    pl.df['ATOM']['b_factor'] = series.values
    print("AFTER\n", pl.df['ATOM']['b_factor'])

    pl.to_pdb(outpath)
    
    return pl.pdb_text


# set groupings

cli.add_command(contact_freqs)
cli.add_command(ligs)
cli.add_command(lig_interactions)
cli.add_command(vis)

vis.add_command(freq)


if __name__ == "__main__":
    cli()
