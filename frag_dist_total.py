# %%
"""
Analyze reaction products across energy-dependent trajectory sets.

This script reads reaction data from simulation folders, maps molecular
identifiers to their corresponding representations, determines the occurrence
percentage of reaction products at each energy, and plots the distributions
of single-, double-, and triple-fragment products.

Expected inputs:
    analysis.pkl
        Pickled dictionary containing molecular identifier information.

    <energy>/<trajectory>/reactions.dat
        Reaction information for each trajectory.

    <energy>/<trajectory>/movie.xyz
        Trajectory file used to determine whether a trajectory should be
        included in the analysis.

Output:
    Frag_dist_total.png
        Figure showing the occurrence percentages of the most common
        fragmentation products as a function of energy.
"""

import numpy as np
import os
import matplotlib.pyplot as plt
from collections import Counter
from rdkit import Chem
import re
import pickle


# Dictionary containing the reaction product associated with each trajectory.
# Structure: reactions[energy][trajectory] = product
reactions = {}

# List containing all observed reaction-product formulas.
all_formulas = []
all_fragment_numbers = set()

# Molecular lookup data loaded from the pickle file.
data = {}


# Load previously generated molecular analysis information.
with open("analysis.pkl", "rb") as f:
    data = pickle.load(f)
    
    # Add an explicit mapping for molecular hydrogen.
    data["[HH]NeighborList"] = ["[H]"]


# Identify directories whose names are integers. These directory names
# represent the different energy sets and are sorted numerically.
set_folders = sorted([folder for folder in os.listdir() if folder.isdigit()], key=int)

# Loop over each energy set.
for set_ in set_folders :
    reactions[set_] = {}

    # Obtain all trajectories belonging to the current energy set.
    trajectories = os.listdir(set_)

    # Process each trajectory independently.
    for trajectory in trajectories :

        # Only process trajectories containing both the trajectory file and
        # the corresponding reaction-analysis file.
        if os.path.isfile(f"{set_}/{trajectory}/reactions.dat") :
            reactions[set_][trajectory] = []

            # Read the reaction information for the current trajectory.
            with open(f"{set_}/{trajectory}/reactions.dat") as f :
                # Use only the final line, corresponding to the final reaction products.
                line = f.readlines()[-1]

                # Extract molecular identifiers from the tab-separated line.
                elements = list(map(str, line.rstrip("\n").split("\t")[1:]))

                # Convert molecular identifiers using the lookup dictionary.
                elements = [data[mol][0] for mol in elements]
                
                all_fragment_numbers.add(len(elements))

                # Store simplified molecular formulas for each product.
                symbols = []

                # Convert each molecular representation into an elemental formula.
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
                        
                    # Count the occurrence of each extracted element.
                    counts = Counter(letters)
                    
                    # Construct an alphabetically ordered molecular formula.
                    symbol = "".join([f"{element}{number}" if number != 1 else f"{element}" for element, number in sorted(counts.items(), key=lambda x : x[0])])
                    symbols.append(symbol)
                    
                # Store single-fragment products directly.
                if len(elements) == 1 :
                    reactions[set_][trajectory] = elements[0]
                        
                    all_formulas.extend(elements)

                # For fragmentation into multiple products, combine the
                # individual formulas using "/" as the separator.
                else :
                    formula = "/".join(sorted(symbols))
                    reactions[set_][trajectory] = formula
                    all_formulas.extend([formula])
                    

"""
Calculate the percentage occurrence of each reaction product at every energy.

The resulting ``points`` dictionary has the structure:

    points[formula][energy] = occurrence_percentage

Each trajectory contributes an equal fraction of the total percentage for its
corresponding energy set.
"""

# Dictionary storing occurrence percentages for every observed product.
points = {}

# Initialize one dictionary for each unique reaction-product formula.
points = {formula : {} for  formula in set(all_formulas)}

# Initialize the occurrence of every product to zero at every energy.
for  formula in all_formulas :
    points[formula] = {}
    for energy, set_ in reactions.items() :  
        points[formula][float(energy)] = 0


# Count the occurrence of each reaction product at every energy.
for energy, set_ in reactions.items() :   
    for trajectory, reaction in set_.items() :  

        # Convert each trajectory occurrence into a percentage of all
        # trajectories belonging to the current energy set.
        if len(reaction.split("/")) == 1 :
            points[reaction][float(energy)] += 1/len(set_)*100
        else :
            points[reaction][float(energy)] += 1/len(set_)*100


"""
Plot the energy-dependent distributions of fragmentation products.

The three panels contain:
    a) The five most frequently observed single-fragment products.
    b) All observed two-fragment products.
    c) All observed three-fragment products.

Chemical labels are formatted using LaTeX-style subscripts, and the final
figure is saved as ``Frag_dist_total.png``.
"""

# Create three vertically stacked panels sharing the same energy axis.
fg, ax = plt.subplots(len(all_fragment_numbers), 1, figsize=(4*len(all_fragment_numbers), 8), sharex=True)

# Select products containing only a single molecular fragment.
single_frag = {species_name : species for species_name, species in points.items() if len(species_name.split("/")) == 1}

# Rank single-fragment products by their total occurrence across all energies
# and retain only the five most frequently observed products.
single_frag = sorted(single_frag.items(), key=lambda x : sum(x[1].values()), reverse=True)[:5]

# Plot each of the five most common single-fragment products.
for formula, species in single_frag :

    # Calculate the total occurrence of the species over all energies.
    counts = sum([i for i in species.values()])

    # Ignore species with zero total occurrence.
    if counts != 0 :
        xpoints = species.keys()
        ypoints = species.values()
        
     
        # Components used to construct a condensed chemical formula.
        condensed_parts = []

        # Interpret the molecular representation as a SMILES string.
        mol = Chem.MolFromSmiles(formula)


        # Ensure implicit hydrogens are calculated
        mol.UpdatePropertyCache()
    
        
        # Loop sequentially through the atoms as they appear in the SMILES
        for atom in mol.GetAtoms():
            symbol = atom.GetSymbol()

            # Determine the total number of hydrogens associated with the atom.
            num_h = atom.GetTotalNumHs()
            
            # Format the individual atom component (e.g., CH3, CH2, CH)
            if num_h == 0:
                h_str = ""
            elif num_h == 1:
                h_str = "H"
            else:
                h_str = f"H{num_h}"

            # Append the condensed representation of the current atom.
            condensed_parts.append(f"{symbol}{h_str}")
        
        # Join all atomic components into a condensed molecular formula.
        e = "".join(condensed_parts)

        # Convert numerical characters into LaTeX-style subscripts.
        label = r"$\mathrm{" + re.sub(r"(\d+)", r"_{\1}", e) + r"}$"

        # Plot the product occurrence as a function of energy.
        ax[0].plot(xpoints, ypoints, "o-", label=label, markersize=6, linewidth=1, alpha=0.8)

# Format the single-fragment panel.
ax[0].tick_params(axis="both", labelsize=16)
ax[0].set_ylim(-0.5, 18)
ax[0].set_ylabel("Occurence [%]", fontsize=16)
ax[0].legend(fontsize=10, loc="upper left")

# Convert energy-folder names to integers for plotting/reference.
xpoints = list(map(int, set_folders))

# Separate products containing two and three fragments.
for i, n_frags in enumerate(sorted(all_fragment_numbers)[1:], 1) :
    current_frags = {species_name : species for species_name, species in points.items() if len(species_name.split("/")) == n_frags}

    # Rank fragmented products by their total occurrence across all energies.
    current_frags = sorted(current_frags.items(), key=lambda x : sum(x[1].values()), reverse=True)
   

    for formula, specie_count in current_frags :

        # Construct and format the product label.
        species_label = "".join(formula)
        label = r"$\mathrm{" + re.sub(r"(\d+)", r"_{\1}", species_label) + r"}$"
        
        # Plot occurrence percentage against energy.
        ax[i].plot(list(specie_count.keys()), list(specie_count.values()), "o-", markersize=6, linewidth=1, alpha=0.8, label=label)



    # Format the two-fragment panel.
    ax[i].tick_params(axis="both", labelsize=16)
    ax[i].set_ylabel("Occurence [%]", fontsize=16)
    ax[i].legend(fontsize=10)
    ax[i].set_title(f"{n_frags} fragments", loc="left", fontsize=20)


# Label the shared energy axis.
ax[-1].set_xlabel("Energy [eV]", fontsize=16)

# Add panel identifiers.
ax[0].set_title("Isomers", loc="left", fontsize=20)

# Adjust subplot spacing to prevent labels and legends from overlapping.
plt.tight_layout()

# Save the completed figure and display it.
fg.savefig("Frag_dist_total.png", bbox_inches="tight")
plt.show()
