# %%
"""
Analyze molecular dynamics trajectories to identify molecular fragments,
determine their chemical connectivity, and calculate fragment-resolved
kinetic properties.

The script reads one or more XYZ trajectory files from a directory specified
on the command line. For every trajectory frame, atomic connectivity is first
represented as a NetworkX graph using an ASE neighbor list with a fixed
distance cutoff. Connected components of this graph are interpreted as
individual molecular fragments.

For each fragment, RDKit is used to determine molecular connectivity and bond
orders and to generate a SMILES representation. If RDKit cannot determine
bond orders, connectivity is instead estimated using ASE's NeighborList and
the resulting structure is represented using single bonds.

Several dynamical properties are calculated for each fragment:

    - Angular momentum
    - Translational kinetic energy
    - Rotational kinetic energy
    - Vibrational/internal kinetic energy

The identified species and fragment properties are written to ``.dat`` files
inside each trajectory directory. Aggregate species populations across all
trajectories are additionally written to the ``species`` directory.

Expected directory structure
----------------------------
The first command-line argument specifies the trajectory-set directory::

    <set>/
        <trajectory_1>/
            movie.xyz
            movie1.xyz
            movie2.xyz
            ...
        <trajectory_2>/
            movie.xyz
            ...

Trajectory files must follow the naming pattern ``movie*.xyz``.

Command-line arguments
----------------------
1. set_folder : str
    Directory containing trajectory subdirectories.

Outputs
-------
For each trajectory directory:

    reactions.dat
        Fragment SMILES identified at every trajectory frame.

    angular_momentum.dat
        Angular momentum magnitude for each fragment.

    translational_energy.dat
        Translational kinetic energy for each fragment.

    rotational_energy.dat
        Rotational kinetic energy for each fragment.

    vibrational_energy.dat
        Vibrational/internal kinetic energy for each fragment.

Additionally:

    species/<set>.dat
        Number of occurrences of each identified species at every frame,
        summed over all trajectories in the set.

Notes
-----
Fragment detection and chemical bond assignment are performed separately.
The initial NetworkX graph determines which atoms belong to the same
fragment, while RDKit subsequently attempts to determine the detailed
chemical bonding within each fragment.
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
from rdkit.Geometry import Point3D
from rdkit.Chem import rdDetermineBonds
from collections import Counter
import re


# The first command-line argument specifies the directory containing
# the trajectory subdirectories to analyze.
set_folders = [sys.argv[1]]



# Dictionary containing all trajectory structures.
#
# Organization:
# structures[set_folder][trajectory] -> list of ASE Atoms frames
structures = {}

# Alternative code for automatically selecting numerically named directories.
#set_folders = sorted([folder for folder in os.listdir() if folder.isdigit()], key=int)

# Loop over all requested trajectory-set directories.
for set_ in set_folders :
    structures[set_] = {}

    # Each item inside the set directory is treated as a potential trajectory.
    trajectories = os.listdir(set_)

    for trajectory in trajectories :

        # Process the trajectory only when its primary movie.xyz file exists.
        if os.path.isfile(f"{set_}/{trajectory}/movie.xyz") :
            atoms = []

            # Find all files matching movie.xyz, movie1.xyz, movie2.xyz, etc.
            # Files are sorted according to the integer appearing in their
            # filename, with movie.xyz treated as the first file.
            movie_files = sorted([file for file in os.listdir(f"{set_}/{trajectory}") if re.search(r"^movie\d*.xyz$", file)], key=lambda x: int(re.search(r"\d+", x).group()) if re.search(r"\d+", x) else 0)

            # Read and concatenate all trajectory segments.
            for i, movie in enumerate(movie_files) :

                # For trajectory segments after the first one, skip frame 0.
                # This avoids including the duplicated boundary frame when
                # consecutive movie files overlap.
                if i > 0 :
                    movie_atoms = read(f"{set_}/{trajectory}/{movie}", index="1:")

                # Keep every frame from the first trajectory file.
                else :
                    movie_atoms = read(f"{set_}/{trajectory}/{movie}", index=":")
                
                # Append the frames from this file to the complete trajectory.
                atoms.extend(movie_atoms)

                print(f"File {set_}/{trajectory}/{movie} was loaded.")

            # Store the complete trajectory.
            structures[set_][trajectory] = atoms


# Dictionary containing connectivity graphs for every trajectory frame.
#
# Organization:
# graphs[set][trajectory][frame] -> NetworkX graph
graphs = {}

for energy, set_ in structures.items() :
    
    graphs[energy] = {}

    for trajectory, structure in set_.items() :
        
        graphs[energy][trajectory] = []

        # Construct an atomic connectivity graph for every trajectory frame.
        for frame in structure :

            # Find pairs of atoms separated by less than the specified
            # 2.5 Å cutoff.
            i, j = neighbor_list('ij', frame, cutoff=2.5, self_interaction=True)

            # Create an undirected graph representing atomic connectivity.
            G = nx.Graph()

            # Add an edge for every neighboring atom pair.
            for atom_i, atom_j in zip(i, j) :
                G.add_edge(atom_i, atom_j)

            # Store the connectivity graph for this trajectory frame.
            graphs[energy][trajectory].append(G)

        print(f"Fragments for trajectory {energy}/{trajectory} have been determined.")


def get_bonds_isomers(graph, atoms) :
    """
    Identify molecular fragments and calculate their dynamical properties.

    Parameters
    ----------
    graph : networkx.Graph
        Atomic connectivity graph for a single trajectory frame. Graph nodes
        correspond to atom indices, and connected components are interpreted
        as separate molecular fragments.

    atoms : ase.Atoms
        Complete atomic configuration corresponding to ``graph``. Atomic
        positions, chemical symbols, momenta, velocities, and masses are
        obtained from this object.

    Returns
    -------
    formulas : list of str
        SMILES representations of all fragments identified in the frame.
        When RDKit bond-order determination fails, the SMILES string is
        generated from ASE NeighborList connectivity using single bonds and
        the suffix ``"NeighborList"`` is appended.

    data : dict
        Fragment-resolved dynamical properties. The dictionary contains:

        ``"angular_momentum"``
            Magnitude of the angular momentum of each fragment.

        ``"translational_energy"``
            Center-of-mass translational kinetic energy of each fragment.

        ``"rotational_energy"``
            Rotational kinetic energy of each fragment.

        ``"vibrational_energy"``
            Remaining kinetic energy after subtracting translational and
            rotational contributions.

    Notes
    -----
    RDKit first attempts to infer molecular connectivity and bond orders
    directly from the three-dimensional atomic coordinates. If bond-order
    determination raises a ``ValueError``, ASE's NeighborList is used as a
    fallback connectivity model.

    For multi-atom fragments, rotational kinetic energy is calculated using
    the angular momentum projected onto the principal axes of inertia.
    Single-atom fragments are assigned zero rotational energy.
    """

    # Each connected component of the graph represents an independent
    # molecular fragment.
    atoms_in_fragments = list(nx.connected_components(graph))

    # List that will contain the SMILES representation of each fragment.
    formulas = []

    # Dictionary containing fragment-resolved dynamical properties.
    data = {}
    data["angular_momentum"] = []
    data["translational_energy"] = []
    data["rotational_energy"] = []
    data["vibrational_energy"] = []
    
    # Store chemical formulas and ASE Atoms objects for all fragments.
    symbols = []
    molecules = []

    # Construct an ASE Atoms object for each connected component.
    for atoms_in_fragment in atoms_in_fragments :
        error = False
        molecule = Atoms()

        # Copy each atom belonging to the current connected component.
        for atom_id in atoms_in_fragment :
            molecule += atoms[atom_id]

        molecules.append(molecule)

        # Store the conventional chemical formula of the fragment.
        symbols.append(molecule.get_chemical_formula())

    # Analyze chemical connectivity and dynamics for every fragment.
    for molecule in molecules :

        # Chemical formula of the current fragment.
        symbol = molecule.get_chemical_formula()

        # By default, attempt RDKit-based bond determination.
        status = True

        # Molecular charge supplied to RDKit.
        charge = 0

        # Attempt to determine molecular connectivity and bond orders
        # directly from the three-dimensional geometry using RDKit.
        if status :

            #fragments.append(molecule)

            # Create an RDKit conformer capable of storing the 3D coordinates
            # of all atoms in the fragment.
            conf = Chem.Conformer(len(molecule))

            # Editable RDKit molecule used to construct the molecular graph.
            rw = Chem.RWMol()

            # Add atoms and their Cartesian coordinates to the RDKit molecule.
            for i, atom in enumerate(molecule) :
                position = atom.position

                conf.SetAtomPosition(i, Point3D(position[0], position[1], position[2]))

                rw.AddAtom(Chem.Atom(atom.symbol))

            # Convert the editable molecule into an RDKit molecule.
            mol = rw.GetMol()

            # Attach the three-dimensional conformer.
            mol.AddConformer(conf, assignId=True)

            try :
                # Create a copy on which RDKit connectivity and bond orders
                # can be assigned.
                conn_mol = Chem.Mol(mol)

                # Determine which atoms are chemically connected based on
                # their positions.
                rdDetermineBonds.DetermineConnectivity(conn_mol, charge=charge) # covFactor=1.6,
                
                # Assign bond orders while allowing charged fragments.
                rdDetermineBonds.DetermineBondOrders(conn_mol, allowChargedFragments=True, charge=charge, embedChiral=False, maxIterations=500)

                # Convert the resulting molecular graph to a SMILES string.
                smiles = Chem.MolToSmiles(conn_mol)
        
            # Mark the RDKit procedure as unsuccessful if bond determination
            # fails.
            except ValueError:
                status = False

        # If RDKit bond determination failed, estimate connectivity using
        # ASE's NeighborList instead.
        if not status :
            print("Switching to NeighborList.")

            # Editable RDKit molecule for the fallback representation.
            rw = Chem.RWMol()

            chemical_symbols = molecule.get_chemical_symbols()

            # Neighbor-list skin distance.
            skin = 0.1

            # Generate element-dependent natural cutoff distances.
            cutoffs = natural_cutoffs(molecule, mult=1.0)

            # Increase the natural cutoffs for hydrogen and carbon.
            for i, atom in enumerate(molecule)  :
                if atom.symbol == "H" :
                    cutoffs[i] *= 1.4

                if atom.symbol == "C" :
                    cutoffs[i] *= 1.1

            # Construct the ASE neighbor list using the modified cutoffs.
            nl_fragment = NeighborList(cutoffs=cutoffs, self_interaction=False, bothways=False, skin=skin)

            nl_fragment.update(molecule)

            # Collect all atom pairs interpreted as bonds.
            bonds = []

            for atom in range(len(molecule)) :

                # Obtain the neighboring atoms of the current atom.
                indices, offsets = nl_fragment.get_neighbors(atom)
                
                for index in indices :
                    bonds.append((atom, index))
            
            # Add all atoms to the RDKit molecule.
            for atom in chemical_symbols :
                rw.AddAtom(Chem.Atom(atom))

            # Represent every neighbor-list connection as a single bond.
            for i, j in bonds :
                rw.AddBond(int(i), int(j), Chem.BondType.SINGLE)

            # Convert the editable molecule into a standard RDKit molecule.
            mol = Chem.Mol(rw)

            # Generate a SMILES representation and mark it as originating
            # from the NeighborList fallback procedure.
            smiles = Chem.MolToSmiles(mol) + "NeighborList"

        # Calculate the magnitude of the fragment's total angular momentum.
        angular_momentum = np.linalg.norm(molecule.get_angular_momentum())

        # Retrieve momenta, velocities, and atomic masses.
        momenta = molecule.get_momenta()
        velocities = molecule.get_velocities()
        masses = molecule.get_masses()

        # Total linear momentum of the fragment.
        P = momenta.sum(axis=0)

        # Center-of-mass translational kinetic energy:
        # E_trans = |P|^2 / (2 M)
        E_trans = P@P/2/masses.sum()

        # Rotational energy is calculated only for fragments containing
        # more than one atom.
        if len(molecule) != 1 :

            # Calculate the center-of-mass velocity.
            v_COM = np.average(velocities, axis=0, weights=masses)

            # Remove center-of-mass translation from the atomic velocities.
            velocities -= v_COM

            # Atomic coordinates relative to the center of mass.
            r = molecule.positions - molecule.get_center_of_mass()

            # Initialize the fragment angular momentum vector.
            J = np.zeros(3)

            # Calculate angular momentum relative to the center of mass:
            # J = sum_i r_i x (m_i v_i)
            for ri, vi, mi in zip(r, velocities, masses):
                J += np.cross(ri, mi*vi)

            # Obtain principal moments of inertia and principal axes.
            I, axes = molecule.get_moments_of_inertia(vectors=True)

            # Project angular momentum onto the principal inertia axes.
            Jp = axes @ J

            # Rotational kinetic energy:
            # E_rot = 1/2 sum_k J_k^2 / I_k
            #
            # Only positive moments of inertia are included.
            E_rot = 0.5*np.sum(Jp[I > 0]**2/I[I > 0])
            

        else :
            # A single atom has no rotational kinetic energy in this
            # decomposition.
            E_rot = 0.0
        
        # Assign the remaining kinetic energy to vibrational/internal motion.
        E_vib = molecule.get_kinetic_energy() - E_trans - E_rot

        # Store the molecular identity and calculated dynamical properties.
        formulas.append(smiles)

        data["angular_momentum"].append(angular_momentum)
        data["translational_energy"].append(E_trans)
        data["rotational_energy"].append(E_rot)
        data["vibrational_energy"].append(E_vib)

    # Return fragment identities and their corresponding properties.
    return formulas, data


# Containers for species identities and fragment-resolved properties.
all_formulas = {}
species = {}
fragment_data = {}

# Organize each dynamical quantity independently.
fragment_data["angular_momentum"] = {}
fragment_data["translational_energy"] = {}
fragment_data["rotational_energy"] = {}
fragment_data["vibrational_energy"] = {}

# Analyze every trajectory set.
for energy, set_ in graphs.items() :   

    species[energy] = {}
    all_formulas[energy] = []

    fragment_data["angular_momentum"][energy] = {}
    fragment_data["translational_energy"][energy] = {}
    fragment_data["rotational_energy"][energy] = {}
    fragment_data["vibrational_energy"][energy] = {}

    # Analyze each trajectory within the current set.
    for trajectory, frames in set_.items() :    

        species[energy][trajectory] = []

        fragment_data["angular_momentum"][energy][trajectory] = []
        fragment_data["translational_energy"][energy][trajectory] = []
        fragment_data["rotational_energy"][energy][trajectory] = []
        fragment_data["vibrational_energy"][energy][trajectory] = []
        
        # Analyze every trajectory frame.
        for i, graph in enumerate(frames) :

            # Determine molecular fragments, their SMILES representations,
            # and their dynamical properties.
            formulas, data = get_bonds_isomers(graph, structures[energy][trajectory][i])
            
            # Store the identified species together with the frame index.
            species[energy][trajectory].append((i, formulas))

            # Add all species from this frame to the global list for this set.
            all_formulas[energy].extend(formulas)
            
            # Store each fragment-resolved dynamical quantity together with
            # the corresponding frame index.
            fragment_data["angular_momentum"][energy][trajectory].append((i, data["angular_momentum"]))

            fragment_data["translational_energy"][energy][trajectory].append((i, data["translational_energy"]))

            fragment_data["rotational_energy"][energy][trajectory].append((i, data["rotational_energy"]))

            fragment_data["vibrational_energy"][energy][trajectory].append((i, data["vibrational_energy"]))

            
        print(f"Species for {energy}/{trajectory} have been determined.")


# Write the molecular species detected in each trajectory to reactions.dat.
for energy, set_ in species.items() :

    for trajectory, frames in set_.items() :

        with open(f"{energy}/{trajectory}/reactions.dat", "w") as f :

            # Each row contains the frame index followed by all fragment
            # SMILES strings detected in that frame.
            for frame in frames :
                f.write(str(frame[0]) + "\t" + "\t".join(frame[1]) + "\n")


# Write each fragment-resolved dynamical property to a separate file.
#
# The variable `data` contains the same property keys generated by
# get_bonds_isomers().
for key in data.keys() :

    for energy, set_ in fragment_data[key].items() :

        for trajectory, frames in set_.items() :

            with open(f"{energy}/{trajectory}/{key}.dat", "w") as f :

                # Each row contains the frame index followed by the values
                # for all fragments present in that frame.
                for frame in frames :
                    
                    f.write(str(frame[0]) + "\t" + "\t".join(map(str, frame[1])) + "\n")


# Create the output directory for aggregated species populations if needed.
if not os.path.isdir(f"species") :
    os.mkdir(f"species") 


# Count each molecular species across all trajectories for every frame.
for k, energy in enumerate(structures.keys()) :

    # Initialize one population time series for every unique species found
    # within the current trajectory set.
    counts = {formula : np.zeros(len(frames)) for formula in set(all_formulas[energy])}

    # Loop over trajectory-frame indices.
    for i in range(len(frames)) :

        # Accumulate species counts over all trajectories.
        for trajectory in species[energy].keys() :

            # Count the number of occurrences of every species in this frame.
            amount = Counter(species[energy][trajectory][i][1])

            # Add the counts from this trajectory to the aggregate population.
            for name, num in amount.items() :

                counts[name][i] += num

    # Optional sorting/limiting of species by total population.
    #counts_ = sorted(counts.items(), key=lambda x: sum(x[1]), reverse=True)[:max_species]
    
    # Use each molecular species as a column in the output array.
    keys = list(counts.keys())

    # Stack species population time series into a two-dimensional integer
    # array: rows correspond to frames and columns correspond to species.
    data = np.column_stack([counts[k] for k in keys]).astype(int)

    # Write the aggregated species populations. The header contains the
    # corresponding SMILES/species identifiers.
    np.savetxt(f"species/{energy}.dat", data, delimiter="\t", header="\t".join(keys), fmt="%d", comments="")
