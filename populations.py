"""
Visualize time-dependent populations of chemical species.

This script reads population data generated for one or more excitation
energies and plots the population of the most abundant chemical species
as a function of simulation time.

For each excitation-energy dataset, species are ranked according to their
total population over the complete trajectory. The selected species are
then plotted using a consistent color across all excitation energies.

The population data are expected to be stored in the ``populations``
directory as tab-separated ``.dat`` files. The first row of each file
contains species names, while the remaining rows contain their populations
at successive simulation steps.

"""

import numpy as np
import matplotlib.pyplot as plt
import os
import argparse


# Configure the command-line interface for the population analysis script.
parser = argparse.ArgumentParser(
        prog="Visualizes time-depdent populations and evalues rate constants of chemical reactions.",
        description="""This script plots and saves time-dependet populations of all species to
                        the <populations> folder. Also it calculates rate constants which 
                        then uses to calcualte populations of given speces that are later compared
                        with the populations from moleculer dynamics simulation."""
    )

# Excitation-energy datasets that should be included in the analysis.
parser.add_argument("-e", "--energies", nargs="+", required=True, type=str, help="Enter the names of the excitation energy folders that will be used as input.")

# Optional limit on the number of highest-population species to plot.
parser.add_argument("-m", "--max", type=int, help="Use this flag if you want set maximum number of trajectories that will be plotted while they are sorted by total population")

# Time conversion factor determined by the simulation timestep and output stride.
parser.add_argument("-t", "--timestep", required=True, type=int, help="Enter product of simulation timestep and stride so one can plot the correct time.")

# Parse command-line arguments supplied by the user.
args = parser.parse_args()

# Maximum number of species displayed for each excitation energy.
max_species = args.max

# Simulation timestep/stride product used for constructing the time axis.
timestep = args.timestep

# Directory containing the population data files.
folder = "populations"

# Excitation-energy datasets requested through the command line.
set_folders = args.energies

# Continue only when the population-data directory exists.
if os.path.isdir(folder):

    # Create one subplot for each excitation-energy dataset.
    # np.atleast_1d ensures that ``ax`` remains indexable even when there
    # is only one excitation-energy dataset.
    fg, ax = plt.subplots(len(set_folders), 1, figsize=(4*len(set_folders), 8), sharex=True)
    ax = np.atleast_1d(ax) 

    # ---------------------------------------------------------
    # First pass: collect the species that will be plotted
    # ---------------------------------------------------------

    # Store unique species selected across all excitation energies.
    species_to_plot = set()

    # Examine each excitation-energy dataset to identify its most
    # populated species.
    for set_ in set_folders:

        # Read species names from the first line of the data file.
        with open(f"{folder}/{set_}.dat", "r") as f:
            header = f.readline().strip().split("\t")

        # Load population values while skipping the species-name header.
        data = np.loadtxt(
            f"{folder}/{set_}.dat",
            delimiter="\t",
            skiprows=1,
            dtype=int
        )

        # Map each species name to its population trajectory.
        # A separate case is required when the file contains only one
        # species because NumPy then returns a one-dimensional array.
        if len(header) == 1:
            d = {header[0]: data}
        else:
            d = {key: data[:, i] for i, key in enumerate(header)}

        # Rank species by their total accumulated population.
        counts_ = sorted(d.items(), key=lambda x: sum(x[1]), reverse=True)

        # Add the highest-population species to the global plotting set.
        species_to_plot.update(name for name, count in counts_[:max_species])

    # ---------------------------------------------------------
    # Assign ONE fixed color to each species
    # ---------------------------------------------------------

    # Retrieve Matplotlib's default plotting color cycle.
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]

    # Assign each species a fixed color so that the same species has
    # the same appearance in every subplot.
    species_colors = {
        species: colors[i % len(colors)]
        for i, species in enumerate(sorted(species_to_plot))
    }

    # ---------------------------------------------------------
    # Plot
    # ---------------------------------------------------------

    # Store a single representative line for each species. These line
    # objects are later used to construct one shared figure legend.
    legend_handles = {}

    # Process and plot each excitation-energy dataset.
    for k, set_ in enumerate(set_folders):

        # Read species names from the header of the current data file.
        with open(f"{folder}/{set_}.dat", "r") as f:
            header = f.readline().strip().split("\t")

        # Load population trajectories for the current excitation energy.
        data = np.loadtxt(
            f"{folder}/{set_}.dat",
            delimiter="\t",
            skiprows=1,
            dtype=int
        )

        # Associate each species with its corresponding population data.
        if len(header) == 1:
            d = {header[0]: data}
        else:
            d = {key: data[:, i] for i, key in enumerate(header)}

        # Create a lookup table mapping each species name to its column index.
        data_ids = {name: i for i, name in enumerate(header)}

        # Rank species by their total population over the trajectory.
        counts_ = sorted(d.items(), key=lambda x: sum(x[1]), reverse=True)

        # Plot only the requested number of highest-population species.
        for name, count in counts_[:max_species]:

            # Convert trajectory indices into simulation time in picoseconds.
            xpoints = np.arange(len(count)) * timestep / 1000
    
            # Plot the population trajectory using the species-specific color.
            line, = ax[k].plot(xpoints, count, label=name, lw=2.5, color=species_colors[name])
            
            # Store only one handle per species for the shared legend.
            if name not in legend_handles:
                legend_handles[name] = line

        # Configure the appearance and labeling of the current subplot.
        ax[k].tick_params(axis="both")
        ax[k].ticklabel_format(axis="y")
        ax[k].set_ylabel("Population")

        # Identify the excitation energy represented by this subplot.
        ax[k].set_title(
            f"Excitation Energy {set_} eV",
            loc="left"
        )

    # ---------------------------------------------------------
    # One legend for the whole figure
    # ---------------------------------------------------------

    # Construct a single shared legend instead of repeating the species
    # legend separately on every subplot.
    leg = fg.legend(
        legend_handles.values(),
        [name for name in legend_handles],
        loc="center left",
        bbox_to_anchor=(0.82, 0.5)
    )

    # Only the bottom subplot requires an x-axis label because all
    # subplots share the same time axis.
    ax[-1].set_xlabel("Time [ps]")


    # Leave space for the shared legend
    # Reserve space on the right for the legend
    fg.tight_layout(rect=[0, 0, 0.80, 1])

    # Display the completed population figure.
    plt.show()
