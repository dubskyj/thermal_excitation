"""
Generate initial conditions for molecular dynamics trajectories.

For each specified excitation energy, the script creates a directory containing
the requested number of trajectory folders. Initial velocities are sampled
along the molecular vibrational normal modes, after which translational and
rotational motion is removed and the momenta are rescaled to match the target
excitation energy.

In restart mode, the script searches existing trajectory directories for movie
files and extracts the final frame of the latest movie file to generate the
next input geometry.
"""

from ase.io import read, write
from fairchem.core import FAIRChemCalculator
from ase.vibrations import Vibrations
import numpy as np
import os
from fairchem.core.units.mlip_unit import load_predict_unit
from ase.md.velocitydistribution import Stationary, ZeroRotation
from ase.units import fs
import argparse
import re


# Define the command-line interface for initial-condition generation.
parser = argparse.ArgumentParser(
        prog="Creates initial conditions for trajectories.",
        description="""This script creates initial conditions for trajectories at given excitation energies.
                    For each excitation energy, a folder is created containing the given number
                    of trajectory folders. In each trajectory folder, a file with initial conditions is created."""
    )

# Specify the molecular geometry used to generate the initial conditions.
parser.add_argument("-i", "--input", required=True, type=str, help="Enter the name of the geometry file that will be used as input.")

# Specify one or more target excitation energies.
parser.add_argument("-e", "--energies", nargs="+", required=True, type=int, help="Enter a sequence of excitation energies for which folders will be created.")

# Specify the number of independent trajectories generated at each excitation energy.
parser.add_argument("-n", "--number", type=int, help="Enter the number of trajectories that will be prepared for each excitation energy.")


parser.add_argument("-m", "--model", required=True, type=str, help="Enter the path to the eSEN model that will be used.")

# Indicate that the input molecule is linear when determining its vibrational degrees of freedom.
parser.add_argument("-l", "--linear", action="store_true", help="Use this flag if the current molecule is linear.")

# Enable restart mode for previously generated trajectories.
parser.add_argument("-r", "--restart", action="store_true", help="Use this flag to create input files for restarting already completed trajectories.")

args = parser.parse_args()


# A trajectory count is required only when generating new initial conditions.
if not args.restart and args.number is None:
    parser.error("--number is required unless --restart is specified.")

# Store command-line arguments in variables used throughout the workflow.
sets = args.energies
n_trajectories = args.number
prepare_for_restart = args.restart
initial_geometry = args.input
output_configuration = "geom.xyz"
is_linear = args.linear
dataset = args.model

# Reduced Planck constant in eV fs.
hbar = 0.6582119569 # eV fs


# Restart mode: generate new input geometries from the latest available trajectory frames.
if prepare_for_restart :
    for set_ in sets :

        # Process the excitation-energy set only if its directory already exists.
        if str(set_) in os.listdir() :

            # Inspect each trajectory belonging to the current excitation-energy set.
            for trajectory in os.listdir(str(set_)) :
                        
                    # Identify movie files and sort them according to their numerical indices.
                    movie_files = sorted([file for file in os.listdir(f"{set_}/{trajectory}") if re.search(r"^movie\d*.xyz$", file)], key=lambda x: int(re.search(r"\d+", x).group()) if re.search(r"\d+", x) else 0)
                    
                    if len(movie_files) != 0 :
                    
                        # Extract the final molecular configuration from the latest movie file.
                        last_frame = read(f"{set_}/{trajectory}/{movie_files[-1]}", index="-1")

                        # Determine the index of the next restart geometry.
                        geom_num = int(re.findall(r"^movie(\d*).xyz$", movie_files[-1])[0]) if re.findall(r"^movie(\d*).xyz$", movie_files[-1])[0] != "" else 0
                        geom_num += 1

                        # Write the extracted configuration as the next trajectory input geometry.
                        write(f"{set_}/{trajectory}/geom{geom_num}.xyz", last_frame)
                        print(f"Input file {set_}/{trajectory}/geom{geom_num}.xyz was generated.")
                    
                    else :
                        print(f"There are no movie files in {set_}/{trajectory}.")
                        
        else :
            print(f"Set {set_} has not been found!")



if not prepare_for_restart :
    # %%

    # Load the machine-learned potential and initialize the ASE calculator.
    predictor = load_predict_unit(dataset)
    calc = FAIRChemCalculator(predictor, task_name="omol")

    # Read the initial molecular geometry and attach the calculator.
    furan = read(initial_geometry)
    furan.calc = calc
    
    # Determine the number of atoms and non-vibrational degrees of freedom.
    # Linear molecules have five non-vibrational modes, whereas nonlinear molecules have six.
    n_atoms = len(furan)
    n_non_vib_nodes = 5 if is_linear else 6
    
    print(f"Number of atoms: {n_atoms}")
    print(f"Vibrational degrees of freedom: {3*n_atoms - n_non_vib_nodes}")
    print(f"Non-vibrational modes: {n_non_vib_nodes})")
    
    # %%
    # Perform the vibrational frequency analysis.
    vib = Vibrations(furan)
    vib.run()
    vib.summary()


    # Obtain the atomic masses for mass weighting of the Hessian.
    masses = furan.get_masses()

    # Construct the inverse square-root mass matrix and diagonalize the mass-weighted Hessian.
    reduced_masses = np.diag(np.repeat(masses**(-1/2), 3))
    eigenvalues, eigenvectors = np.linalg.eigh(reduced_masses @ vib.H @ reduced_masses)

    # Keep only vibrational modes.
    eigenvalues = eigenvalues[n_non_vib_nodes:]

    # Convert the eigenvalues to angular frequencies in fs^-1.
    omega = np.sqrt(eigenvalues) * np.sqrt(1/103.642696)

    # Transform the normal modes into Cartesian coordinates.
    modes = eigenvectors.T.reshape(3*n_atoms, n_atoms, 3)[n_non_vib_nodes:]
    cart_modes = np.diag(masses**(-1/2)) @ modes


    # Generate an independent set of trajectories for each excitation energy.
    for set in sets :
        print("-"*40)

        # Create the excitation-energy directory or continue numbering trajectories
        # from the highest existing trajectory index.
        if not os.path.isdir(f"{set}") :
            os.mkdir(f"{set}") 
            max_id = 0
            print(f"Excitation energy folder {set} has been created.")
        else :
            max_id = int(sorted(os.listdir(f"{set}"), key=int)[-1]) + 1
            print(f"Excitation energy folder {set} already exists.")
        
        print("-"*40)
        os.path.isdir(f"{set}")

        # Generate the requested number of independent trajectories.
        for trajectory in np.arange(max_id, n_trajectories + max_id) :

            # Create a directory for the current trajectory if it does not already exist.
            if not os.path.isdir(f"{set}/{trajectory}") :
                os.mkdir(f"{set}/{trajectory}") 

            # Set the target kinetic energy to the specified excitation energy.
            E_target = set

            # Calculate the effective mass associated with each Cartesian normal mode.
            mu = np.array([1/np.sum(mode**2) for mode in cart_modes])

            # Initialize normal-mode momenta and Cartesian atomic velocities.
            p = np.zeros(len(omega))
            vel = np.zeros((len(masses), 3))

            # Sample the momentum of each vibrational mode from its zero-point
            # Gaussian distribution and accumulate its contribution to the atomic velocities.
            for i in range(len(omega)):
                sigma_q = np.sqrt(hbar/(2*mu[i]*omega[i]))
                sigma_p = np.sqrt(hbar*mu[i]*omega[i]/2)

                p[i] = np.random.normal(0.0, sigma_p)

                vel  += p[i]/np.sqrt(mu[i])*cart_modes[i]

            # Assign the sampled velocities to the molecular geometry.
            furan.set_velocities(vel/fs)
            #print(furan.get_kinetic_energy(), np.sum(p**2 / (2*mu)))

            # Remove net translational and rotational motion.
            Stationary(furan)
            ZeroRotation(furan)

            # Rescale the momenta so that the kinetic energy matches the target excitation energy.
            Ek = furan.get_kinetic_energy()
            furan.set_momenta(furan.get_momenta() * np.sqrt(E_target/Ek)) 

            # Write the prepared initial configuration to the trajectory directory.
            write(f"{set}/{trajectory}/{output_configuration}", furan)

            print(f"Input file {set}/{trajectory}/{output_configuration} was generated.")

    print("#### Done ####")
