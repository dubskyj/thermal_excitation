"""
Analyze the conservation of total energy and total angular momentum in
molecular dynamics trajectories.

The script reads one or more sets of molecular dynamics trajectories stored
in XYZ files. Each input directory is interpreted as a set of trajectories
associated with a particular excitation energy.

For every trajectory, the script:

    1. Locates trajectory files named ``movie.xyz``, ``movie1.xyz``,
       ``movie2.xyz``, etc.
    2. Loads and concatenates the trajectory frames.
    3. Calculates the total energy for every frame.
    4. Calculates the magnitude of the total angular momentum for every frame.
    5. Computes the total-energy drift relative to the initial frame.
    6. Plots the energy drift and angular momentum as functions of simulation
       time.

When a trajectory is split across multiple movie files, the first frame of
each subsequent file is skipped. This assumes that it duplicates the final
frame of the preceding file.

The effective timestep supplied through ``--timestep`` corresponds to the
time interval between consecutive stored trajectory frames in femtoseconds.
The plotted time axis is converted from femtoseconds to picoseconds.


Expected directory structure
----------------------------
energy_directory/
    trajectory_1/
        movie.xyz
        movie1.xyz
        movie2.xyz
        ...
    trajectory_2/
        movie.xyz
        movie1.xyz
        ...

"""

import numpy as np
import os
from ase.io import read
import matplotlib.pyplot as plt
import re
import argparse


# Configure the command-line interface.
parser = argparse.ArgumentParser(prog="Analysis of total energy and total angular momentum conservation.", description="Analyze and visualize the conservation of total energy and total angular momentum along molecular dynamics trajectories.")

parser.add_argument("-e", "--energies", nargs="+", required=True, type=str, help="Enter the names of the excitation energy folders that will be used as input.")
parser.add_argument("-t", "--timestep", required=True, type=float, help="Enter the effective time interval between consecutive trajectory frames in fs (simulation timestep × output stride).")

# Parse command-line arguments.
args = parser.parse_args()

# Store the input trajectory directories and effective timestep.
set_folders = args.energies
timestep = args.timestep


data = {}



# =============================================================================
# Load trajectory data
# =============================================================================

# Loop over all supplied energy directories.
for set_ in set_folders:
    data[set_] = {}

    # Obtain and sort all entries in the current energy directory.
    trajectories = sorted(os.listdir(set_))

    # Process each trajectory independently.
    for trajectory in trajectories:
        trajectory_path = os.path.join(set_, trajectory)

        # Ignore entries that are not directories.
        if not os.path.isdir(trajectory_path):
            continue

        # Locate trajectory files following the naming convention
        # movie.xyz, movie1.xyz, movie2.xyz, ...
        #
        # The numerical component is used to ensure that split trajectory
        # files are processed in chronological order.
        movie_files = sorted([file for file in os.listdir(trajectory_path) if re.search(r"^movie\d*\.xyz$", file)], key=lambda x: int(re.search(r"\d+", x).group()) if re.search(r"\d+", x) else 0)

        # Skip directories that do not contain any trajectory files.
        if len(movie_files) == 0:
            continue

        # Container for all ASE Atoms objects belonging to this trajectory.
        frames = []

        # Read and concatenate all trajectory segments.
        for i, movie in enumerate(movie_files):
            movie_path = os.path.join(trajectory_path, movie)

            # Read every frame from the first trajectory file.
            if i == 0:
                movie_atoms = read(movie_path, index=":")
            # For subsequent files, skip the first frame because it is
            # assumed to duplicate the final frame of the previous file.
            else:
                movie_atoms = read(movie_path, index="1:")

            # Append the newly loaded frames to the complete trajectory.
            frames.extend(movie_atoms)

            # Report successful loading of each trajectory segment.
            print(f"File {movie_path} was loaded.")

        # Store the calculated total energy and magnitude of the total angular
        # momentum for every frame in the trajectory.
        data[set_][trajectory] = {}
        data[set_][trajectory]["total_energy"] = np.array([frame.get_total_energy() for frame in frames])
        data[set_][trajectory]["total_angular_momentum"] = np.array([np.linalg.norm(frame.get_angular_momentum()) for frame in frames])


# %%
# =============================================================================
# Plot colors
# =============================================================================

# Use the tab20 categorical color map to distinguish trajectories.
# Colors are reused if an energy set contains more than 20 trajectories.
colors = plt.cm.tab20.colors



# %%
# =============================================================================
# Plot all trajectories together for each excitation energy
# =============================================================================

# Create one row of plots for each supplied energy directory.
#
# Left column:  total-energy drift relative to the initial frame.
# Right column: magnitude of the total angular momentum.
#
# squeeze=False guarantees that ax remains a two-dimensional array even when
# only one energy directory is supplied.
fg, ax = plt.subplots(len(set_folders), 2, figsize=(12, 4*len(set_folders)), sharex=True, squeeze=False)

# Loop over the energy sets stored in the data dictionary.
for i, (energy, trajectories) in enumerate(data.items()):
    print(i, energy)


    # Plot every trajectory belonging to the current energy set.
    for j, (trajectory, data_) in enumerate(trajectories.items()):
        total_energy = data_["total_energy"]
        total_angular_momentum = data_["total_angular_momentum"]

        # Determine the trajectory length directly from the number of
        # available energy values.
        n_steps = len(total_energy)

        # Construct the simulation-time axis and convert fs to ps.
        time = np.arange(n_steps) * timestep / 1000

        # Calculate the energy change relative to the initial trajectory frame.
        energy_drift = total_energy - total_energy[0]

        # Assign a different color to each trajectory. The modulo operation
        # cycles through the available colors when necessary.
        color = colors[j % len(colors)]

        # Plot the total-energy drift and angular-momentum magnitude.
        ax[i][0].plot(time, energy_drift, "o", markersize=0.4, color=color)
        ax[i][1].plot(time, total_angular_momentum, "o", markersize=0.4, color=color)

    # Configure the total-energy subplot for the current excitation energy.
    ax[i][0].set_title(f"    Excitation Energy {energy} eV", loc="left")
    ax[i][0].set_ylabel(r"$E(t)-E(0)$ [eV]", fontsize=12)
    ax[i][0].tick_params(axis="both", labelsize=10)
    ax[i][0].ticklabel_format(axis="y", style="sci", scilimits=(0, 0))


    # Configure the angular-momentum subplot for the current excitation energy.
    ax[i][1].set_ylabel(r"$|\mathbf{L}|$", fontsize=12)
    ax[i][1].tick_params(axis="both", labelsize=10)
    ax[i][1].ticklabel_format(axis="y", style="sci", scilimits=(0, 0))

# Add time-axis labels only to the bottom row because the x-axis is shared
# between all subplot rows.
ax[-1][0].set_xlabel("Time [ps]", fontsize=12)
ax[-1][1].set_xlabel("Time [ps]", fontsize=12)

# Automatically adjust subplot spacing and display the completed figure.
plt.tight_layout()
plt.show()
