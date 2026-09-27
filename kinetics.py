"""
Construct a kinetic model from reaction trajectories and compare the resulting
species populations with populations obtained directly from the trajectories.

The script loads molecular identifiers, processes reaction trajectories,
smooths short-lived changes in molecular identity, identifies reaction-state
transitions, calculates transition rates, constructs a rate matrix, propagates
the kinetic model in time, and converts configuration populations into
individual species populations.

The resulting kinetic information and reaction rates are saved for subsequent
analysis.
"""

import numpy as np
import matplotlib.pyplot as plt
import os
from collections import Counter
from scipy.linalg import expm
import pickle
import argparse


"""
Load reaction trajectories for the selected excitation-energy set.

Reaction identifiers stored in ``reactions.dat`` are converted using the
molecular mapping contained in ``analysis.pkl``. The resulting ``reactions``
dictionary contains the molecular species present at every recorded frame of
each trajectory.
"""

parser = argparse.ArgumentParser(
        prog="Visalizes time-depdent populations and evalues rate constants of chemical reactions.",
        description="""This script plots and saves time-dependet populations of all species to
                        the <populations> folder. Also it calculates rate constants which 
                        then uses to calcualte populations of given speces that are later compared
                        with the populations from moleculer dynamics simulation."""
    )

parser.add_argument("-e", "--energy", required=True, type=str, help="Enter the name of the geometry file that will be used as input.")
parser.add_argument("-m", "--max", type=int, help="Use this flag if you want set maximum number of trajectories that will be plotted while they are sorted by total population")
parser.add_argument("-r", "--rates", action="store_true", help="Use this flag if you want to calculate rate constants and see the comparsion between molecular dynamics trajecotires and calculated trajecotires.")
parser.add_argument("-s", "--steps", required=True, type=int, help="Enter the number of steps that is common for all the trajecotires.")
parser.add_argument("-t", "--timestep", required=True, type=float, help="Enter product of simulation timestep and stride so one can plot the correct time.")

args = parser.parse_args()


# Time interval associated with consecutive stored trajectory frames.
timestep = args.timestep

excitation_energy = args.energy

# Maximum number of species intended for plotting.
max_species = args.max

calculate_rate_constants = args.rates

# Directory containing the direct species population data.
folder = "populations"

# Dictionary containing molecular species for every trajectory and frame.
reactions = {}

# Excitation-energy folders included in the analysis.
set_folders = [excitation_energy] 


# Load the molecular identifier mapping generated during previous analysis.
with open("analysis.pkl", "rb") as f:
    data = pickle.load(f)
    
    # Add an explicit mapping for molecular hydrogen.
    data["[HH]NeighborList"] = ["[H]"]


for set_ in set_folders :
    reactions[set_] = {}

    trajectories = os.listdir(set_)
    for trajectory in trajectories :

        if os.path.isfile(f"{set_}/{trajectory}/reactions.dat") :

            with open(f"{set_}/{trajectory}/reactions.dat") as f :
                lines = f.readlines().copy()
                
                if len(lines) > args.steps :
                    reactions[set_][trajectory] = []
                    
                    for line in lines :
                        elements = list(map(str, line.rstrip("\n").split("\t")[1:]))

                        reactions[set_][trajectory].append(elements)
                        #reactions[set_][trajectory].append([data[mol][0] for mol in elements])



            #print(len(reactions[set_][trajectory]))
            #print(f"File {set_}/{trajectory}/movie.xyz was loaded.")

# Create a compact state name for every unique molecular representation.
names = {}

# Collect all unique molecular representations contained in the lookup data.
all_names = set([i[0] for i in data.values()])

# Assign sequential state identifiers S0, S1, S2, ...
for i, line in enumerate(sorted(all_names)) :
    names[line] = "S" + str(i)

"""
Smooth the molecular-state trajectories.

A state change is accepted only when the species detected in the current frame
remain consistent over a number of subsequent frames defined by
``steps_ahead``. Temporary or inconsistent changes are replaced by the
previous accepted state.

States containing ``.`` are also rejected because this notation indicates
that RDKit has interpreted the molecular representation as multiple fragments.
"""

# Smooths the species trajectory

# Number of future frames used to determine whether a state change is stable.
steps_ahead = 4

# Dictionary containing the smoothed trajectories.
trajs = {}

# Process trajectories separately for every excitation energy.
for energy, set_ in reactions.items() :

    trajs[energy] = {}

    # Smooth each trajectory independently.
    for trajectory, reaction in set_.items() :

        # Copy the original sequence of molecular states.
        frames = [a for a in reaction]

        # Initialize the smoothed trajectory using the first frame.
        traj = [frames[0]]

        # Examine every subsequent frame.
        for i in range(1, len(frames)) :
            if i < len(frames) :

                # Accept the current state when the molecular species present
                # remain consistent within the forward-looking frame window.
                if set(list(set(np.concatenate(frames[i:i+steps_ahead])))) == set(np.array(frames[i])) :
                    traj.append(frames[i])
                else :

                    # Otherwise retain the previously accepted state.
                    traj.append(traj[i-1])

                # Every specie with "." which means that RDkit defines that as framgnets
                if any("." in element for element in traj[i]) :

                    # Reject states interpreted by RDKit as disconnected fragments.
                    traj[i] = traj[i-1]

        # Store the smoothed trajectory.
        trajs[energy][trajectory] = traj


"""
Convert each smoothed trajectory into a sequence of state-change events.

Consecutive frames containing the same collection of molecular species are
collapsed into a single state. Each stored entry therefore contains the frame
index at which a new molecular configuration first appears and the species
present in that configuration.
"""

# Shortens the blocks of same species
sequences = {}

# Collect every molecular formula appearing in the processed trajectories.
all_formulas = set()

# Process every smoothed trajectory.
for energy, set_ in trajs.items() :

    sequences[energy] = {}

    for trajectory, traj in set_.items() :

        frames = traj

        # The first state always begins at frame zero.
        products = [(0, frames[0])]

        # Add a new state whenever the molecular composition changes.
        for i in range(1, len(frames)) :
            if frames[i-1] != frames[i] :
                products.append((i, frames[i]))
    
    
        """
        Optional filtering of temporary fragmentation events.

        This disabled block identifies intermediate states containing more
        fragments than both the preceding and following states. Such events
        correspond to patterns of the form A + B -> C + D + E -> A + B and
        can be removed from the state sequence.

        to_remove = []
        for j, product in enumerate(products[1:-1], 1) :
            prev, current, next = products[j-1][1], products[j][1], products[j+1][1]

            # Filtres out the reactions such A + B -> C + D + E -> A + B
            if len(current) > 1 :
                
                if len(prev) == len(next) and len(prev) < len(current) and len(next) < len(current) :
                    to_remove.append(j)

        for k in sorted(to_remove, reverse=True):
            products.pop(k)
        """

        # Collect every molecular formula encountered in the state sequence.
        for product in products :
            all_formulas.update(product[1])    

        # Store the compressed sequence of state changes.
        sequences[energy][trajectory] = products



"""
Calculate and save direct species populations from the smoothed trajectories.

For every frame, the number of occurrences of each molecular species across
all trajectories is counted. The resulting population matrix is written to a
tab-separated file for the corresponding excitation energy.
"""

# Generate the species population files when the output directory does not exist.
for k, energy in enumerate(set_folders) :
    
    if not os.path.isdir(folder) :
        os.mkdir(folder)

    # Initialize a time-dependent count array for every observed species.
    counts = {formula : np.zeros(len(frames)) for formula in set(all_formulas)}

    # Count each species across all trajectories at every frame.
    for i in range(len(frames)) :
        for trajectory in trajs[energy].keys() :

            # Count the molecular species in the current trajectory frame.
            amount = Counter(trajs[energy][trajectory][i])

            # Add these counts to the total species populations.
            for name, num in amount.items() :
                counts[name][i] += num

    
    # Convert the species population dictionary into a matrix.
    keys = list(counts.keys())
    data = np.column_stack([counts[k] for k in keys]).astype(int)

    # Save species names as the header and populations as integer columns.
    np.savetxt(f"{folder}/{energy}.dat", data, delimiter="\t", header="\t".join(keys), fmt="%d", comments="")
     


"""
Load and plot the directly counted species populations.

Population files generated from the smoothed trajectories are read and plotted
as a function of simulation time.
"""

# Plot species populations when the population directory exists.
if os.path.isdir(folder) :

    # Create the population figure.
    fg, ax = plt.subplots(1, 1, figsize=(10, 3*len(set_folders) + 1), sharex=True)

    # Process each excitation-energy dataset.
    for k, set_ in enumerate(set_folders) :

        # Read molecular species names from the file header.
        with open(f"{folder}/{set_}.dat", "r") as f:
            header = f.readline().strip().split("\t")

        # Load the corresponding population matrix.
        data = np.loadtxt(f"{folder}/{set_}.dat", delimiter="\t", skiprows=1, dtype=int)

        # Convert the loaded data into a species-indexed dictionary.
        if len(header) == 1 :
            d = {header[0]: data}
        else :
            d = {key: data[:, i] for i, key in enumerate(header)}
            

        # Rank species according to their integrated population.
        counts_ = sorted(d.items(), key=lambda x: sum(x[1]), reverse=True)

        # Plot every species population as a function of time.
        for name, count in counts_[:max_species] :
            xpoints = np.arange(len(count))*timestep
            ax.plot(xpoints, count, label=name)

        # Format the population plot.
        #ax.legend()
        ax.set_title(f"{set_} eV", fontsize=12)
        ax.set_ylabel("Count", fontsize=12)

    # Label the shared time axis.
    ax.set_xlabel("Time [ps]", fontsize=12)



# Display the direct population plot.
leg = fg.legend(
    fontsize=8,
    loc="center left",
    bbox_to_anchor=(0.82, 0.5)
)

fg.tight_layout(rect=[0, 0, 0.80, 1])
plt.show()





if calculate_rate_constants :

    """
    Determine state residence times and observed transitions.

    For each trajectory, the time spent in every molecular configuration is
    accumulated. Whenever one configuration changes into another, a textual
    transition of the form

        reactants -> products

    is stored for subsequent calculation of transition rates.
    """

    # Total residence time associated with each molecular configuration.
    time = {}

    # Observed transitions for every excitation energy.
    transitions = {}

    for energy, set_ in sequences.items():

        time[energy] = {}
        transitions[energy] = []

        # Process each compressed trajectory sequence.
        for trajectory, sequence in set_.items():

            # Examine each molecular state in the sequence.
            for i, frame in enumerate(sequence):

                # Starting frame of the current state.
                t0 = frame[0]

                # Construct a reproducible configuration name by sorting species.
                state = " + ".join(sorted(frame[1]))

                # For the final state, use the final trajectory frame as its endpoint.
                if i + 1 == len(sequence):
                    t1 = len(reactions[energy][trajectory]) - 1   # or times[-1] if using real time
                else:

                    # Otherwise the state ends when the next configuration begins.
                    t1 = sequence[i + 1][0]

                    # Construct and record the observed state-to-state transition.
                    products = " + ".join(sorted(sequence[i + 1][1]))
                    transition = f"{state} -> {products}"
                    transitions[energy].append(transition)
                
                # Accumulate the residence time of the current state.
                time[energy][state] = time[energy].get(state, 0) + (t1 - t0)

        
     
    """
    Construct the kinetic transition-rate matrix.

    Observed transition counts are divided by the total residence time of their
    reactant configurations and by ``timestep``. These rates populate the
    off-diagonal elements of the generator matrix ``Q``. The diagonal elements
    are subsequently set to the negative sum of all outgoing rates.
    """

    # Number of unique molecular configurations in the selected energy set.
    N = len(time[excitation_energy])

    # Initialize the kinetic generator matrix.
    Q = np.zeros((N, N))

    # Store calculated transition rates for each excitation energy.
    rates = {}

    for energy, tranistions_set in transitions.items() :


        rates[energy] = {}

        # Count how many times each transition occurs.
        counts = Counter(tranistions_set)
      

        # Assign each molecular configuration an integer matrix index.
        ids = {key : i for i, key in enumerate(time[energy].keys())}

        # Calculate a rate for every observed transition.
        for reaction, count in counts.items() :

            # Separate the initial and final configurations.
            reactants, products = reaction.split(" -> ")

            # Obtain their corresponding matrix indices.
            i = ids[reactants]
            j = ids[products]

            # Calculate the transition rate from event count and residence time.
            rate = count/time[energy][reactants]/timestep

            # Populate the off-diagonal rate matrix and store the value.
            Q[i, j] = rate
            rates[energy][reaction] = rate

    # Set each diagonal element to minus the total outgoing rate from that state.
    for i in range(len(Q)):
        Q[i, i] = -np.sum(Q[i, :])



    # %%

    # Visualize the kinetic transition-rate matrix.
    plt.title("Transition-rate  matrix")
    plt.imshow(Q)
    plt.show()


    # %%
    """
    Construct the stoichiometric mapping between molecular configurations and
    individual species.

    Each column of ``S`` represents one molecular configuration, while each row
    represents one molecular species. Matrix elements contain the number of times
    a species occurs in the corresponding configuration.
    """

    # List all molecular configurations included in the kinetic model.
    configs = list(time[energy].keys())

    # Map each configuration to an integer index.
    cfg_ids = {cfg : i for i, cfg in enumerate(configs)}

    # Extract every unique molecular species occurring in the configurations.
    species = sorted({
        molecule
        for config in configs
        for molecule in config.split(" + ")
    })


    # Assign each molecular species an integer index.
    ids_ = {s:i for i, s in enumerate(species)}

    # Initialize the species-to-configuration stoichiometric matrix.
    S = np.zeros((len(species), len(configs)))

    # Populate the stoichiometric matrix.
    for j, cfg in enumerate(configs) :

        # Separate the individual molecular species in the configuration.
        cfg_ = cfg.split(" + ")

        # Count each species appearing in the configuration.
        for sp in cfg_ :
            i = ids_[sp]
            S[i][j] += 1

    # Indicate completion of the stoichiometric matrix.
    print("S matrix done!")

    # Visualize the species/configuration mapping.
    plt.title("Stoichiometric matrix")
    plt.imshow(S.T)
    plt.show()

    # %%
    """
    Propagate configuration populations using the kinetic rate matrix.

    The initial population is placed entirely in the specified parent molecular
    configuration. The matrix exponential of the transposed generator matrix is
    then used to obtain configuration populations at each requested time.
    """

    # Use the direct-population time points for kinetic propagation.
    times = xpoints

    # Initialize the configuration population vector.
    P0 = np.zeros(N)

    # Place all initial population in the parent molecular configuration.
    P0[ids["[H]C1=C([H])C([H])=C([H])O1"]] = len(reactions[excitation_energy]) 

    # Allocate storage for configuration populations at every time point.
    Pcfg = np.zeros((len(times), len(configs)))

    # Propagate the kinetic model to each requested time.
    for i,t in enumerate(times):

        # Evaluate the continuous-time Markov-model population.
        Pcfg[i] = expm(Q.T*t) @ P0
        
        print(i, t)

    
    # ==========================================================
    # Convert into species populations
    # ==========================================================

    # Convert configuration populations into individual molecular-species populations.
    Pspecies = (S @ Pcfg.T).T
    


    # %%
    """
    Compare direct trajectory populations with populations predicted by the
    kinetic model and save the modeled populations for later analysis.
    """

    # Dictionary used to retain the directly measured population data.
    popul = {}

    # Load and plot direct trajectory populations when available.
    if os.path.isdir(folder) :
        popul[folder] = {}

        # Create a figure for the directly measured species populations.
        fg, ax = plt.subplots(1, 1, figsize=(10, 3*len(set_folders) + 1), sharex=True)

        # Process each selected excitation energy.
        for k, set_ in enumerate(set_folders) :
            
            
            # Read species names from the population-file header.
            with open(f"{folder}/{set_}.dat", "r") as f:
                header = f.readline().strip().split("\t")

            # Load the direct population data.
            data = np.loadtxt(f"{folder}/{set_}.dat", delimiter="\t", skiprows=1, dtype=int)

            # Retain the population matrix for later use.
            popul[folder][set_] = data
            
            # Convert the population matrix into a species-indexed dictionary.
            if len(header) == 1 :
                d = {header[0]: data}
            else :
                d = {key: data[:, i] for i, key in enumerate(header)}

            # Rank species according to their integrated population.
            counts_ = sorted(d.items(), key=lambda x: sum(x[1]), reverse=True)

            # Plot each directly measured species population.
            for name, count in counts_ :
                xpoints = np.arange(len(count))*timestep
                ax.plot(xpoints, count, label=name)
            
                
            # Format the direct population plot.
            #ax.legend()
            ax.set_title(f"{set_} eV", fontsize=12)
            ax.set_ylabel("Count", fontsize=12)

        # Label the time axis.
        ax.set_xlabel("Time [ps]", fontsize=12)


        # Plot the species populations predicted by the kinetic model.
        for i,s in enumerate(species) :
            plt.plot(
                times,
                Pspecies[:,i],
                lw=2,
                label=s,
            )
            
        # Format and display the modeled species-population plot.
        # Display the direct population plot.
        leg = fg.legend(
            fontsize=8,
            loc="center left",
            bbox_to_anchor=(0.82, 0.5)
        )

        fg.tight_layout(rect=[0, 0, 0.80, 1])
        plt.xlabel("Time")
        plt.ylabel("Species population")
        plt.tight_layout()
        plt.show()


    # Save the modeled species populations and species-index mapping.
    #with open(f"data_kinetics_{set_folders[0]}.pkl", "wb") as f:
    #    pickle.dump([Pspecies, ids_], f)

    # Create the output directory for transition rates when necessary.
    if not os.path.isdir("rates") :
        os.mkdir("rates")

    # Save the calculated reaction rates for each excitation energy.
    for folder in set_folders :

        with open(f"rates/{folder}.dat", "w") as f :

            # Convert full configuration transitions into net molecular reactions.
            for reaction, rate in rates[folder].items() :

                # Separate the initial and final configurations.
                state_i, state_j = reaction.split(" -> ")

                # Convert each configuration into lists of molecular species.
                left = state_i.split(" + ")
                right = state_j.split(" + ")

                # Determine which species are consumed and produced.
                lost = Counter(left) - Counter(right)
                gained = Counter(right) - Counter(left)

                # Construct the net reactant and product strings.
                reactants = " + ".join(lost)
                products = " + ".join(gained)

                # Write and print the resulting net reaction and transition rate.
                f.write(f"{reactants} -> {products}\tk={rate*1000} ps^-1\n")
                print(f"{reactants} -> {products}\tk={rate*1000} ps^-1")
