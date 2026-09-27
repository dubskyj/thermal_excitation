"""
Identify molecular fragments and chemical species from molecular dynamics trajectories.

The script loads trajectory files from excitation-energy directories and constructs
a graph representation of every trajectory frame based on interatomic distances.
Connected components of each graph are interpreted as individual molecular fragments.

For each fragment, RDKit is used to determine molecular connectivity and bond orders
from its Cartesian geometry and to generate a SMILES representation. If RDKit cannot
determine the bonding pattern directly, an ASE neighbor-list approach is used as a
fallback to construct the molecular connectivity.

Geometries corresponding to each identified chemical species are collected across
all excitation energies and trajectories and written to separate XYZ files.
"""

import numpy as np
import os
from ase.io import read
from ase import Atoms
import networkx as nx
from ase.neighborlist import NeighborList, natural_cutoffs
from ase.neighborlist import neighbor_list
from rdkit import Chem
import sys
from rdkit.Chem import rdDetermineBonds
from collections import Counter
from rdkit.Geometry import Point3D
from ase.io import write
import re

# %%
# Store molecular dynamics trajectories according to excitation energy and trajectory.
structures = {}

# Identify excitation-energy directories and sort them numerically.
set_folders = sorted([folder for folder in os.listdir() if folder.isdigit()], key=int)

# Load all available trajectory frames from each excitation-energy directory.
for set_ in set_folders :
    structures[set_] = {}

    trajectories = os.listdir(set_)
    for trajectory in trajectories :

       # Process only trajectory directories containing the initial movie file.
       if os.path.isfile(f"{set_}/{trajectory}/movie.xyz") :
            atoms = []

            # Identify all movie files and sort them according to their numerical indices.
            movie_files = sorted([file for file in os.listdir(f"{set_}/{trajectory}") if re.search(r"^movie\d*.xyz$", file)], key=lambda x: int(re.search(r"\d+", x).group()) if re.search(r"\d+", x) else 0)

            # Load consecutive movie files while avoiding duplication of the first frame
            # in continuation files.
            for i, movie in enumerate(movie_files) :
                if i > 0 :
                    movie_atoms = read(f"{set_}/{trajectory}/{movie}", index="1:")
                else :
                    movie_atoms = read(f"{set_}/{trajectory}/{movie}", index=":")
                
                atoms.extend(movie_atoms)
                print(f"File {set_}/{trajectory}/{movie} was loaded.")

            # Store the complete sequence of trajectory frames.
            structures[set_][trajectory] = atoms


# %%
# Construct graph representations of all molecular dynamics frames.
graphs = {}

for energy, set_ in structures.items() :
    
    graphs[energy] = {}
    for trajectory, structure in set_.items() :
        
        graphs[energy][trajectory] = []

        # Represent each frame as a graph in which atoms are nodes and neighboring
        # atoms within the specified cutoff are connected by edges.
        for frame in structure :
            i, j = neighbor_list('ij', frame, cutoff=2.5, self_interaction=True)
            
            G = nx.Graph()
            for atom_i, atom_j in zip(i, j) :
                G.add_edge(atom_i, atom_j)

            graphs[energy][trajectory].append(G)

        print(f"Fragments for trajectory {energy}/{trajectory} have been determined.")


# %%
def get_bonds_isomers(graph, atoms) :
    """
    Identify molecular fragments and determine their chemical connectivity.

    Connected components of the molecular graph are interpreted as individual
    fragments. For each fragment, its Cartesian geometry is transferred to an
    RDKit molecule and RDKit is used to determine connectivity and bond orders.
    The resulting molecular species is represented by a SMILES string.

    If RDKit cannot determine a valid bonding pattern from the geometry, an ASE
    neighbor list is used to construct the connectivity and all detected bonds
    are represented as single bonds.

    Parameters
    ----------
    graph
        NetworkX graph describing the connectivity of atoms in the complete
        molecular system.
    atoms
        ASE Atoms object containing the atomic species and Cartesian coordinates
        corresponding to the graph.

    Returns
    -------
    formulas
        List of SMILES representations for the molecular fragments.
    geoms
        List of ASE Atoms objects containing the corresponding fragment geometries.
    """

    # Determine individual molecular fragments from the connected components
    # of the molecular graph.
    atoms_in_fragments = list(nx.connected_components(graph))
    formulas = []
    geoms = []

    
    symbols = []
    molecules = []

    # Construct a separate ASE Atoms object for every connected fragment.
    for atoms_in_fragment in atoms_in_fragments :
        error = False
        molecule = Atoms()
        for atom_id in atoms_in_fragment :
            molecule += atoms[atom_id]

        molecules.append(molecule)
        symbols.append(molecule.get_chemical_formula())

    # Determine the chemical connectivity and SMILES representation of each fragment.
    for molecule in molecules :
        symbol = molecule.get_chemical_formula()
        status = True
        charge = 0


        if status :
            #fragments.append(molecule)

            # Transfer the Cartesian coordinates and atomic identities to an
            # RDKit molecular representation.
            conf = Chem.Conformer(len(molecule))

            rw = Chem.RWMol()
            for i, atom in enumerate(molecule) :
                position = atom.position
                conf.SetAtomPosition(i, Point3D(position[0], position[1], position[2]))
                rw.AddAtom(Chem.Atom(atom.symbol))

            mol = rw.GetMol()
            mol.AddConformer(conf, assignId=True)

            # Attempt to determine molecular connectivity and bond orders directly
            # from the three-dimensional geometry using RDKit.
            try :
                conn_mol = Chem.Mol(mol)
                rdDetermineBonds.DetermineConnectivity(conn_mol, charge=charge) # covFactor=1.6,


                rdDetermineBonds.DetermineBondOrders(conn_mol, allowChargedFragments=True, charge=charge, embedChiral=False, maxIterations=500)

                smiles = Chem.MolToSmiles(conn_mol)
        
            # Mark unsuccessful RDKit bond-order assignments for fallback processing.
            except ValueError:
                status = False

        # If RDKit bond determination fails, construct connectivity using an
        # ASE neighbor list.
        if not status :
            print("Switching to NeighborList.")
            rw = Chem.RWMol()
            chemical_symbols = molecule.get_chemical_symbols()

            # Generate element-dependent natural cutoff radii for bond detection.
            skin = 0.1
            cutoffs = natural_cutoffs(molecule, mult=1.0)

            # Increase hydrogen and carbon cutoffs to improve connectivity detection.
            for i, atom in enumerate(molecule)  :
                if atom.symbol == "H" :
                    cutoffs[i] *= 1.4
                if atom.symbol == "C" :
                    cutoffs[i] *= 1.1

            # Construct the neighbor list used to identify bonded atom pairs.
            nl_fragment = NeighborList(cutoffs=cutoffs, self_interaction=False, bothways=False, skin=skin)
            nl_fragment.update(molecule)

            # Collect all bonded atom pairs detected by the neighbor list.
            bonds = []
            for atom in range(len(molecule)) :
                indices, offsets = nl_fragment.get_neighbors(atom)
                
                for index in indices :
                    bonds.append((atom, index))
            
            # Create the corresponding RDKit atoms.
            for atom in chemical_symbols :
                rw.AddAtom(Chem.Atom(atom))

            # Represent all neighbor-list connections as single bonds.
            for i, j in bonds :
                rw.AddBond(int(i), int(j), Chem.BondType.SINGLE)

            # Generate a SMILES representation and mark species identified using
            # the neighbor-list fallback.
            mol = Chem.Mol(rw)
            smiles = Chem.MolToSmiles(mol) + "NeighborList"


        # Store the molecular identifier and its corresponding geometry.
        formulas.append(smiles)
        geoms.append(molecule)


    return formulas, geoms


# %%
# Collect identified molecular species and their corresponding geometries.
all_formulas = {}
species = {}
geoms = {}


# Analyze every frame of every trajectory across all excitation energies.
for energy, set_ in graphs.items() :   
    for trajectory, frames in set_.items() :    

        for i, graph in enumerate(frames) :

            # Identify all molecular fragments present in the current frame.
            formulas, data = get_bonds_isomers(graph, structures[energy][trajectory][i])

            # Group geometries according to their molecular SMILES representation.
            for formula, geom in zip(formulas, data) :
                if geoms.get(formula) is not None :
                    geoms[formula].append(geom)
                else :
                    geoms[formula] = [geom]
                

# %%
# Create the output directory for geometries if it does not already exist.
if not os.path.isdir(f"geoms") :
    os.mkdir(f"geoms") 


# Write all geometries associated with each identified chemical species
# to a separate XYZ file.
for formula, data in geoms.items() :
    write(f"geoms/{formula}.xyz", data)

    

