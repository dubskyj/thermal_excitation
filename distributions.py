"""
Analyze the energy and angular-momentum distributions of reaction products.

The script reads reaction products and fragment-resolved observables from
trajectory directories at different excitation energies. Reaction products
are converted into molecular formulas, and the final translational,
rotational, vibrational, and angular-momentum values are grouped according
to fragment identity.

The resulting distributions are visualized as stacked histograms for each
excitation energy. A final scatter plot compares translational and rotational
energies of the individual fragments.

Expected trajectory files:
    movie.xyz
        Molecular dynamics trajectory.

    reactions.dat
        Reaction-product information.

    angular_momentum.dat
        Fragment angular momenta.

    translational_energy.dat
        Fragment translational energies.

    rotational_energy.dat
        Fragment rotational energies.

    vibrational_energy.dat
        Fragment vibrational energies.

Outputs:
    Energy_distributions.png
        Stacked distributions of translational, rotational, and vibrational
        energies.

    Ang_distributions.png
        Stacked angular-momentum distributions.
"""


import numpy as np
import os
import matplotlib.pyplot as plt
import re
from collections import Counter
from matplotlib.patches import Patch
import sys


"""
Read reaction products and fragment-resolved observables from all trajectory
directories.

The ``reactions`` dictionary stores the molecular fragments produced by each
trajectory, while ``data`` stores the corresponding angular momentum and
energy observables. ``all_formulas`` records the unique fragment formulas
observed at each excitation energy.
"""

# Fragment-resolved observables to read from the trajectory directories.
observables = ["angular_momentum", "translational_energy", "rotational_energy", "vibrational_energy"]

# Reaction products indexed by excitation energy and trajectory.
reactions = {}

# Unique molecular formulas observed at each excitation energy.
all_formulas = {}

# Dictionary containing all fragment-resolved observables.
data = {}
data["angular_momentum"] = {}
data["translational_energy"] = {}
data["rotational_energy"] = {}
data["vibrational_energy"] = {}

# Identify numerically named folders and sort them according to excitation energy.
set_folders = sys.argv[1:] #sorted([folder for folder in os.listdir() if folder.isdigit()], key=int)

# Loop over each excitation-energy directory.
for set_ in set_folders :

    # Temporary list containing all fragment formulas found at this energy.
    all_formulas_ = []

    # Initialize dictionaries for the current excitation energy.
    reactions[set_] = {}
    data["angular_momentum"][set_] = {}
    data["translational_energy"][set_] = {}
    data["rotational_energy"][set_] = {}
    data["vibrational_energy"][set_] = {}

    # Obtain all trajectory directories belonging to the current energy.
    trajectories = os.listdir(set_)

    # Process each trajectory independently.
    for trajectory in trajectories :

        # Only include trajectories containing a molecular-dynamics trajectory file.
        if os.path.isfile(f"{set_}/{trajectory}/movie.xyz") :


            # Initialize reaction and observable storage for this trajectory.
            reactions[set_][trajectory] = []
            data["angular_momentum"][set_][trajectory] = []
            data["translational_energy"][set_][trajectory] = []
            data["rotational_energy"][set_][trajectory] = []
            data["vibrational_energy"][set_][trajectory] = []

            # Read reaction-product information when available.
            if os.path.isfile(f"{set_}/{trajectory}/reactions.dat") :
                with open(f"{set_}/{trajectory}/reactions.dat") as f :

                    # Use the final line of the reaction file, corresponding to
                    # the final products of the trajectory.
                    line = f.readlines()[-1] #for line in f.readlines() :
                        
                    # Extract the molecular representations from the tab-separated line.
                    elements = list(map(str, line.rstrip("\n").split("\t")[1:]))
        
                    # Store simplified molecular formulas for each fragment.
                    symbols = []

                    # Convert every molecular representation into an elemental formula.
                    for element in elements :

                        # Remove the NeighborList suffix when present.
                        if "NeighborList" in element :
                            element = element.rstrip("NeighborList")
                            
                        # Treat molecular hydrogen as a special case.
                        if element == "[HH]" :
                            letters = ["H"]
                        else :

                            # Extract alphabetic characters representing atomic symbols.
                            letters = [letter for letter in element if letter.isalpha()]
                            
                        # Count occurrences of each extracted element.
                        counts = Counter(letters)
                        
                        # Construct an alphabetically ordered molecular formula.
                        symbol = "".join([f"{element}{number}" if number != 1 else f"{element}" for element, number in sorted(counts.items(), key=lambda x : x[0])])
                        symbols.append(symbol)
                        
                    # Store the fragment formulas associated with this trajectory.
                    reactions[set_][trajectory] = symbols

                    # Add the formulas to the complete list for this excitation energy.
                    all_formulas_.extend(symbols)
                    
                   
                    
                    # Alternative handling of reaction products retained from
                    # the original implementation.
                    #if len(elements) != 1 :
                    #    reactions[set_][trajectory] = elements
                    #    all_formulas_.extend(elements)

            # Read each requested fragment-resolved observable.
            for key in observables :

                # Only process the observable if its corresponding file exists.
                if os.path.isfile(f"{set_}/{trajectory}/{key}.dat") :
                    with open(f"{set_}/{trajectory}/{key}.dat") as f :

                        # Read every recorded timestep from the observable file.
                        for line in f.readlines() :
    
                            # Ignore the first column and convert the remaining
                            # fragment values to floating-point numbers.
                            elements = list(map(float, line.rstrip("\n").split("\t")[1:]))
                    
                            # Store the observable values for this timestep.
                            data[key][set_][trajectory].append(elements)             

    # Store only the unique fragment formulas found at this excitation energy.
    all_formulas[set_] = set(all_formulas_)
    



# %%
"""
Group final fragment properties according to molecular formula.

For every excitation energy and trajectory, the final recorded angular
momentum, translational energy, rotational energy, and vibrational energy
of each fragment are assigned to the corresponding molecular formula.

The resulting structure is:

    points[observable][energy][formula] = [values]
"""

# Dictionary containing fragment-resolved distributions for each observable.
points = {}
points["angular_momentum"] = {}
points["translational_energy"] = {}
points["rotational_energy"] = {}
points["vibrational_energy"] = {}

# Loop over each excitation energy.
for energy, set_ in reactions.items() :   

    # Initialize one list per molecular formula for each observable.
    points["angular_momentum"][energy] = {formula : [] for  formula in all_formulas[energy]}
    points["translational_energy"][energy] = {formula : [] for  formula in all_formulas[energy]}
    points["rotational_energy"][energy] = {formula : [] for  formula in all_formulas[energy]}
    points["vibrational_energy"][energy] = {formula : [] for  formula in all_formulas[energy]}

    # Process each trajectory belonging to the current excitation energy.
    for trajectory, reaction in set_.items() :    
        

        # Assign each fragment's final observables to its molecular formula.
        for i, formula in enumerate(reaction) :
            points["angular_momentum"][energy][formula].append(data["angular_momentum"][energy][trajectory][-1][i])
            points["translational_energy"][energy][formula].append(data["translational_energy"][energy][trajectory][-1][i])
            points["rotational_energy"][energy][formula].append(data["rotational_energy"][energy][trajectory][-1][i])
            points["vibrational_energy"][energy][formula].append(data["vibrational_energy"][energy][trajectory][-1][i])


# Create a grid containing one row per selected excitation energy and one
# column for translational, rotational, and vibrational energy.
fg, ax = plt.subplots(len(reactions), 3, figsize=(12, 8), sharex="col", squeeze=False)    

# Human-readable axis labels for each energy observable.
y_labels = {"translational_energy" : "Translational Energy [eV]",
            "rotational_energy" : "Rotational Energy [eV]",
            "vibrational_energy" : "Vibrational Energy [eV]"}

# Maximum histogram energies used to define the plotting ranges.
max_energy = {"translational_energy" : 2.6,
            "rotational_energy" : 2,
            "vibrational_energy" : 5.5}


# Get all species appearing in the plots
# Hydrogen and the intact C4H4O species are excluded from the plotted products.
all_species = {
    formula
    for key in ["translational_energy", "rotational_energy",
                "vibrational_energy", "angular_momentum"]
    for energy in points[key]
    for formula, values in points[key][energy].items()
    if len(values) != 0 and formula not in ["H", "C4H4O"]
}


colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]

species_colors = {
    formula: colors[i % len(colors)]
    for i, formula in enumerate(all_species)
}

# Plot translational, rotational, and vibrational energy distributions.
for i, key in enumerate(["translational_energy", "rotational_energy", "vibrational_energy"]) :
    
    # Skip the first three energy entries and plot the remaining excitation energies.
    for j, (energy, species_) in enumerate(list(points[key].items())) :

        # Rank species according to the number of available data points.
        species = sorted(species_.items(), key=lambda x: len(x[1]), reverse=True)
        
        # Obtain colors for species that contain data and are included in the plot.
        colors = [
        species_colors[formula]
        for formula, data_points in species
        if len(data_points) != 0 and formula not in ["H", "C4H4O"]
            ]

        
        
        # Collect the individual species distributions and their labels.
        datapoints = []
        label = []

        for k, (formula, data_points) in enumerate(species) :

            # Exclude empty distributions, atomic hydrogen, and intact C4H4O.
            if len(data_points) != 0 and formula != "H" and formula != "C4H4O" :

                # Convert the values to a NumPy array for histogram plotting.
                datapoints.append(np.array(data_points))

                # Format numerical stoichiometric coefficients as LaTeX subscripts.
                formula_updated = r"$\mathrm{" + re.sub(r"(\d+)", r"_{\1}", formula) + r"}$" 
                label.append(formula_updated)


        # Plot normalized stacked histograms for all fragment species.
        ax[j][i].hist(datapoints, bins=np.linspace(0, max_energy[key], 40), label=label, density=True, histtype="bar", stacked=True, color=colors)

        # Format subplot axes.
        #ax[j][i].set_ylabel(y_labels[key], fontsize=20)
        ax[j][i].tick_params(axis="both", labelsize=25)
        #ax[j][i].legend(fontsize=26)
        #ax[j][i].legend(loc='center left', bbox_to_anchor=(1, 0.5))
      
        # Add panel identifiers and excitation-energy labels to the first row.
        ax[j][1].set_title(f"Excitation Energy {energy} eV", loc="left", fontsize=20)

        # Label the histogram intensity axis for each row.
        ax[j][0].set_ylabel("Intensity", fontsize=20)

# Label each energy-distribution column.
ax[-1][0].set_xlabel("Translational Energy [eV]", fontsize=20)
ax[-1][1].set_xlabel("Rotational Energy [eV]", fontsize=20)
ax[-1][2].set_xlabel("Vibrational Energy [eV]", fontsize=20)

# Construct common legend entries using the predefined species colors.
legend_handles = [
    Patch(facecolor=species_colors[formula], label=r"$\mathrm{" + re.sub(r"(\d+)", r"_{\1}", formula) + r"}$" )
    for formula in all_species
]

fg.legend(
    handles=legend_handles,
    loc="center left",
    bbox_to_anchor=(0.82, 0.5),
    ncol=1,
    fontsize=16
)


# Adjust subplot spacing and save the energy-distribution figure.
plt.tight_layout(rect=[0, 0, 0.8, 1])
fg.savefig("Energy_distributions.png", bbox_inches="tight")
plt.show()

# Create vertically stacked panels for angular-momentum distributions.
fg, ax = plt.subplots(len(reactions), 1, figsize=(12, 8), sharex=True)    

# Select angular momentum as the observable to plot.
key = "angular_momentum"

# Plot the angular-momentum distribution for each selected excitation energy.
for j, (energy, species_) in enumerate(list(points[key].items())) :

    # Rank species according to the number of available data points.
    species = sorted(species_.items(), key=lambda x: len(x[1]), reverse=True)
    
    # Retrieve the predefined colors for the species included in the plot.
    colors = [
        species_colors[formula]
        for formula, data_points in species
        if len(data_points) != 0 and formula not in ["H", "C4H4O"]
    ]

    
    """
    Alternative angular-momentum visualization using individual unfilled
    step histograms. This section is retained but disabled.

    for formula, data_points in species :
        if len(data_points) != 0 and formula != "H" and formula != "C4H4O" :
            ax[j].hist(data_points, bins=40, histtype="step", density=True, fill=False, label=formula)
   """
            
    # Collect angular-momentum distributions and formatted species labels.
    datapoints = []
    label = []

    for k, (formula, data_points) in enumerate(species) :

        # Exclude empty distributions, atomic hydrogen, and intact C4H4O.
        if len(data_points) != 0 and formula != "H"  and formula != "C4H4O" :

            # Convert the values to a NumPy array for histogram plotting.
            datapoints.append(np.array(data_points))

            # Format numerical stoichiometric coefficients as LaTeX subscripts.
            formula_updated = r"$\mathrm{" + re.sub(r"(\d+)", r"_{\1}", formula) + r"}$" 
            label.append(formula_updated)


    # Plot the normalized stacked angular-momentum histogram.
    ax[j].hist(datapoints, bins=np.linspace(0, 10.5, 40), label=label, density=True, histtype="bar", stacked=True, color=colors)          
        
            
    # Format the current angular-momentum subplot.
    ax[j].set_ylabel("Intensity", fontsize=20)
    ax[j].tick_params(axis="both", labelsize=25)
    #ax[j].legend(fontsize=25)

    # Add panel identifiers and excitation-energy labels.
    ax[j].set_title(f"   Excitation Energy {energy} eV", loc="left", fontsize=20)


# Label the shared angular-momentum axis.
ax[-1].set_xlabel(r"Angular Momentum [amu$\frac{\text{Å}^2}{\mathrm{fs}}$]", fontsize=20)

# Construct legend entries using the same species-color mapping as the
# energy-distribution figure.
legend_handles = [
    Patch(facecolor=species_colors[formula], label=r"$\mathrm{" + re.sub(r"(\d+)", r"_{\1}", formula) + r"}$" )
    for formula in all_species
]

# Place the common species legend outside the angular-momentum panels.
fg.legend(
    handles=legend_handles,
    loc="center left",
    bbox_to_anchor=(0.82, 0.5),
    ncol=1,
    fontsize=16
)


# Adjust subplot spacing and save the angular-momentum figure.
plt.tight_layout(rect=[0, 0, 0.8, 1])
fg.savefig("Ang_distributions.png", bbox_inches="tight")
plt.show()

