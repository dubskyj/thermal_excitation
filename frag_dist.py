"""
Analyze molecular dynamics trajectories at different excitation energies.

The script reads trajectory specific reaction data from folders named by
excitation energy, classifies the final molecular state of each trajectory,
and calculates the percentage of trajectories corresponding to different
reaction outcomes and numbers of fragments.

The results are visualized as two stacked
bar charts showing:

1. The relative occurrence of intact molecules, isomerization,
   fragmentation, and fragmentation accompanied by H/H2 loss.
2. The relative occurrence of one, two, or three molecular fragments.

Expected directory structure
----------------------------
<energy>/
    <trajectory>/
        movie.xyz
        reactions.dat

where ``<energy>`` is an integer-valued folder name.
"""

import numpy as np
import os
import matplotlib.pyplot as plt
from collections import Counter

# Store the number of fragments and reaction status for each trajectory.
# The first dictionary level corresponds to excitation energy and the
# second level to individual trajectories.
frags = {}
statuses = {}


# Identify folders whose names are integers (excitation energies) and
# sort them numerically
set_folders = sorted(
    [folder for folder in os.listdir() if folder.isdigit()],
    key=int
)

n_fragments = set()

# Process one excitation-energy folder at a time.
for set_ in set_folders:
    frags[set_] = {}
    statuses[set_] = {}

    # Each subdirectory is assumed to correspond to one trajectory.
    trajectories = os.listdir(set_)

    for trajectory in trajectories:

        # Only process trajectories for which both the molecular
        # trajectory and reaction-analysis files are available.
        if os.path.isfile(f"{set_}/{trajectory}/reactions.dat"):

            with open(f"{set_}/{trajectory}/reactions.dat") as f:
                lines = f.readlines()
                
                
                # The final line describes the species present at the
                # end of the trajectory.
                line = lines[-1]
                
                first_species = lines[0].rstrip("\n").lstrip("0\t")

                # Ignore the first tab-separated field and retain the
                # molecular species reported in the remaining fields.
                elements = list(map(str, line.rstrip("\n").split("\t")[1:]))

                symbols = []

                # Convert each molecular representation into a simplified
                # elemental-composition string. This is used below to
                # identify trajectories involving hydrogen loss.
                for element in elements :

                    # Remove a NeighborList marker if one is present.
                    if "NeighborList" in element:
                        element = element.rstrip("NeighborList")

                    # Extract alphabetic characters from the molecular
                    # representation and count their occurrences.
                    letters = [letter for letter in element if letter.isalpha()]
                    counts = Counter(letters)

                    # Construct an alphabetically ordered composition
                    # string, omitting the number when the count is one.
                    symbol = "".join(
                        [
                            f"{element}{number}"
                            if number != 1
                            else f"{element}"
                            for element, number in sorted(
                                counts.items(),
                                key=lambda x : x[0]
                            )
                        ]
                    )
                    symbols.append(symbol)
                    

                n_fragments.add(len(elements))
                
                # More than one final species indicates fragmentation.
                if len(elements) != 1:

                    # Distinguish fragmentation accompanied by H2/H loss
                    # from other fragmentation pathways.
                    if "H2" in symbols:
                        status = "H_loss"
                    else:
                        status = "frag"

                # A single final species corresponds either to the
                # original molecule or to an isomerized structure.
                else:

                    # Reference representation of the intact molecule.
                    if elements[0] == first_species :
                        status = "intact"
                    else:
                        status = "isom"

                # Record the classification and number of final fragments
                # for this trajectory.
                statuses[set_][trajectory] = status
                frags[set_][trajectory] = len(elements)


# Calculate the percentage of trajectories belonging to each reaction
# category at every excitation energy.
bars = {
    "intact": [],
    "isom": [],
    "frag": [],
    "H_loss": []
}


for energy, set_ in statuses.items():

    # Count trajectories in each reaction category at this energy.
    counts = Counter(set_.values())

    for i in bars.keys():

        # Convert trajectory counts to percentages.
        if counts.get(i):
            bars[i].append(counts[i] / len(set_) * 100)
        else:
            bars[i].append(0)


bars_frags = {i : [] for i in n_fragments}

for energy, set_ in frags.items() :

    # Count trajectories according to their number of final fragments.
    counts = Counter(set_.values())

    for i in bars_frags.keys() :
        if counts.get(i) :
            bars_frags[i].append(counts[i]/len(set_)*100)
        else :
            bars_frags[i].append(0)


# Human-readable labels for the reaction categories used in the plot.
alias = {
    "intact": "Intact",
    "isom": "Isomerization",
    "frag": "Fragmentation",
    "H_loss": r"Fragmentation + H$_2$/H loss"
}

# Generate singular/plural labels for the number of fragments.
alias_frag = lambda num: (
    f"{num} Fragments" if num > 1 else f"{num} Fragment"
)

# Width of each stacked bar.
width = 0.8

# Convert excitation-energy folder names to integers for the x-axis.
xpoints = list(map(int, set_folders))


# Create two vertically stacked panels sharing the excitation-energy axis.
fig, ax = plt.subplots(
    2,
    1,
    figsize=(12, 8),
    sharex=True
)


# -------------------------------------------------------------------------
# Panel (a): reaction outcome
# -------------------------------------------------------------------------

# Track the cumulative height required for stacked bars.
bottom = np.zeros(len(set_folders))

for num, count in bars.items():

    ax[0].bar(
        xpoints,
        count,
        width,
        label=alias[num],
        bottom=bottom
    )

    # The next category is stacked on top of the current categories.
    bottom += count

ax[0].set_ylabel("Ocurrence [%]", fontsize=32)
ax[0].tick_params(axis="both", labelsize=25)
ax[0].legend(fontsize=16)


# -------------------------------------------------------------------------
# Panel (b): number of fragments
# -------------------------------------------------------------------------

# Reset the cumulative height for the second stacked bar chart.
bottom_frags = np.zeros(len(set_folders))

for num_frags, count_frags in bars_frags.items():

    ax[1].bar(
        xpoints,
        count_frags,
        width,
        label=alias_frag(num_frags),
        bottom=bottom_frags
    )

    bottom_frags += count_frags

ax[1].set_ylabel("Ocurrance [%]", fontsize=32)
ax[1].tick_params(axis="both", labelsize=25)
ax[1].legend(fontsize=16)


# Add panel labels.
ax[0].set_title("Main reaction channels", loc="left", fontsize=20)
ax[1].set_title("Number of fragmetation events", loc="left", fontsize=20)

# The excitation-energy axis is shared by both panels.
ax[1].set_xlabel("Excitation Energy [eV]", fontsize=32)

# Optimize spacing and display the figure.
plt.tight_layout()
plt.show()

