"""
Identify representative molecular structures (medoids) from geometry datasets.

This script processes molecular geometry files stored in the ``geoms``
directory and selects a representative structure from each dataset based on
pairwise root-mean-square deviation (RMSD).

For each geometry dataset, the script:

1. Determines the number of structures contained in the XYZ file.
2. Reduces the number of structures considered for particularly large
   datasets.
3. Calculates pairwise RMSD values using optimal rotational alignment with
   the Kabsch algorithm.
4. Computes the mean RMSD of each structure relative to the other structures.
5. Identifies the structure with the lowest mean RMSD as the medoid.
6. Writes the selected medoid to the ``medoids`` directory as an XYZ file.

The RMSD calculation removes translational differences by centering the
structures and rotational differences by applying the optimal Kabsch
rotation. Reflections are explicitly prevented.

Input
-----
geoms/
    Directory containing molecular geometry files in XYZ format.

Output
------
medoids/
    Directory containing the representative medoid structure for each
    geometry dataset.

Notes
-----
The script assumes that structures being compared contain the same number
of atoms and that corresponding atoms occur in the same order.
"""



import numpy as np
import os
from ase.io import iread
from ase.io import write
import subprocess



def rmsd(atoms1, atoms2):
    """
    Calculate the root-mean-square deviation (RMSD) between two structures.

    The function compares the Cartesian coordinates of two ASE `Atoms`
    objects after translating both structures to their respective centers
    of mass-like coordinate centers (the arithmetic mean of all atomic
    positions). A Kabsch alignment is then used to determine the optimal
    rotation between the two structures while preventing reflections.

    Parameters
    ----------
    atoms1 : ase.Atoms
        First atomic structure to compare.
    atoms2 : ase.Atoms
        Second atomic structure to compare.

    Returns
    -------
    float
        The RMSD between the two aligned structures in the same distance
        units as the input coordinates, typically Angstroms.

    Raises
    ------
    ValueError
        If the two structures contain different numbers of atoms.
    """

    # Copy the Cartesian coordinates so that the original ASE structures
    # are not modified during the RMSD calculation.
    P = atoms1.get_positions().copy()
    Q = atoms2.get_positions().copy()

    # Check same number of atoms
    if len(P) != len(Q):
        raise ValueError("Structures have different numbers of atoms")

    # Center both molecules by subtracting the average position of all atoms.
    # This removes any influence of the absolute translation of the structures.
    P -= P.mean(axis=0)
    Q -= Q.mean(axis=0)

    # Kabsch optimal rotation:
    # Construct the covariance matrix between the two centered structures
    # and use singular value decomposition (SVD) to determine the rotation
    # that minimizes the squared coordinate differences.
    H = P.T @ Q
    U, S, Vt = np.linalg.svd(H)

    # Prevent reflection:
    # The Kabsch procedure can sometimes produce an improper rotation
    # (a reflection). Flip the final singular vector when necessary so
    # that the resulting transformation has determinant +1.
    if np.linalg.det(Vt.T @ U.T) < 0:
        Vt[-1, :] *= -1

    # Construct the optimal rotation matrix from the SVD components.
    R = Vt.T @ U.T

    # Rotate the first structure into the coordinate frame of the second.
    P_aligned = P @ R.T

    # Calculate the RMSD as the square root of the mean squared distance
    # between corresponding atoms in the aligned structures.
    return np.sqrt(np.mean(np.sum((P_aligned - Q)**2, axis=1)))



# Create the output directory for medoid structures if it does not already exist.
if not os.path.isdir(f"medoids") :
    os.mkdir(f"medoids") 

# Iterate over every geometry dataset stored in the "geoms" directory.
for i, dataset in enumerate(os.listdir("geoms")) :

    # Remove the final four characters of the filename (e.g. ".xyz")
    # to obtain the molecular formula used for the output filename.
    formula = dataset[:-4]
    print(f" {formula} ".center(60, "#"))
   

    # Count the total number of lines in the geometry file.
    # The value is later used to determine how many structures are present.
    result = subprocess.run(
        ["wc", "-l", f"geoms/{dataset}"],
        capture_output=True,
        text=True,
        check=True
    )
    
    n_lines = int(result.stdout.split()[0])
    
    
    # Read the number of atoms from the first line of the geometry file.
    with open(f"geoms/{dataset}", "r") as f:
        n_atoms = int(f.readline())

    # Determine the number of structures in the file.
    # Each XYZ structure is assumed to contain n_atoms + 2 lines:
    # one atom-count line, one comment line, and one line per atom.
    n = n_lines/(n_atoms+2)
        
    
    
    
    print(f"Structures: {n}")
    
    print(f"Waiting for {dataset}.")

    # Reduce the number of structures loaded for very large datasets.
    # Larger datasets are sampled to reduce memory and computation requirements.
    if n > 2e5 :
        structures = iread(f"geoms/{dataset}", index="::20")
    elif n > 5e4 :
        structures = iread(f"geoms/{dataset}", index="::1")
    else :
        structures = iread(f"geoms/{dataset}", index=":")
 
    print(f"{dataset} loaded!")
 
    # Count the number of structures available after the selected reduction.
    # The iterator is consumed by this operation.
    n = sum(1 for _ in structures)
    print(f"Structures after reduction: {n}")
    
    # Allocate a square matrix to store the pairwise RMSD values
    # between all structures in the reduced dataset.
    rmsd_matrix = np.zeros((n, n))

    print(f"Calculating RMSD.")

    # Compare every structure against every other structure.
    # The resulting RMSD values are stored symmetrically in the matrix.
    for s1, i in enumerate(structures):
        for s2, j in enumerate(structures):
            value = rmsd(s1, s2)

            rmsd_matrix[i, j] = value
            rmsd_matrix[j, i] = value

    # Calculate the mean RMSD of each structure relative to all structures.
    # A lower value indicates a structure that is, on average, more
    # representative of the dataset.
    mean_rmsd = rmsd_matrix.mean(axis=1)

    # Select the structure with the smallest mean RMSD.
    # This structure is the medoid of the dataset under the RMSD metric.
    medoid_idx = np.argmin(mean_rmsd)

    # Read the selected medoid structure from the original geometry file.
    medoid =  iread(f"geoms/{dataset}", index=medoid_idx) #structures[medoid_idx]

    #medoids[formula] = medoid

    # Save the selected medoid as an XYZ file in the "medoids" directory.
    write(f"medoids/{formula}.xyz", medoid)
    print(f"Medoid has been obtained!")

    # Report the mean RMSD of the selected medoid to all structures.
    print("Mean RMSD:", mean_rmsd[medoid_idx], "Å")

