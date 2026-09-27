# Summary

This set of scripts enables the calculation and analysis of trajectories of thermally (vibrationally) excited molecules with arbitrary excitation energies specified in eV.

```text
    Input Geometry
        ↓
    get_sets.py
    Initial Conditions
        ↓
    omol.py / r.omol
    MD Trajectories
        ↓
    reactions.py / r.reactions
    Trajectory Analysis
        ↓
    sort.py
        ↓
    medoids.py
        ↓
    opt.py
        ↓
    analysis.pkl
        ↓
    ├── thermal.py ──────────────> Energy & Angular-Momentum Conservation
    ├── frag_dist.py ────────────> Fragmentation Statistics
    ├── frag_dist_total.py ──────> Detailed Fragmentation Analysis
    ├── distributions.py ────────> Energy & Angular-Momentum Distributions
    ├── plot_reactions.py ───> Reaction Networks
    └── kinetics.py ─────────> Populations & Rate Constants
                ↓
            populations.py
                ↓
        Population Comparison
```

# Initiation and Simulation

First, the initial conditions for the simulations need to be prepared. This is done using the `get_sets.py` script, which creates `<energy>` folders containing the specified number of `<trajectory>` subfolders. Each trajectory folder contains a `geom.xyz` file with the initial atomic positions and velocities.

The `get_sets.py` script accepts command-line arguments. A guide to the available options can be displayed by running:

```bash
py get_sets.py -h
```

A typical command is, for example:

```bash
py get_sets.py -i furan.xyz -e 8 9 10 -n 10
```

In this example, the script reads the molecular geometry from `furan.xyz` and creates folders corresponding to excitation energies of 8, 9, and 10 eV. Each energy folder contains 10 trajectories.

For each trajectory, a file containing the initial conditions is prepared. This file contains the positions of the optimized molecule and the corresponding momenta. The specified excitation energy represents the excess energy relative to the ground-state energy of the optimized structure. This energy is distributed among all vibrational modes, while the momenta are sampled from the Wigner distribution and subsequently scaled so that the total kinetic energy matches the desired value. Finally, the algorithm removes any overall rotational and translational motion, ensuring that the initial state contains only vibrational excitation.

Next, the simulations need to be run. Trajectory integration is performed using the `omol.py` script, which is submitted to the compute nodes using the `r.omol` submission script. The latter contains all required input parameters. The input file corresponds to the `geom.xyz` file described above, while the output file should be set to `movie.xyz`, which contains all saved frames along the trajectory.

In addition, the integration timestep, total number of simulation steps, and output interval (stride) must be specified. The stride determines how frequently trajectory frames are saved. Importantly, the path to the dataset used by the current machine-learned interatomic potential must also be provided.

Once the simulations are complete, further analyses are performed using the `reaction.py` script. These analyses determine the chemical species encountered along the trajectory (`reactions.dat`), the partitioning of energy (`rotational_energy.dat`, `translational_energy.dat`, and `vibrational_energy.dat`), and the angular momentum (`angular_momentum.dat`). The `reaction.py` script is submitted to a computational cluster using the `r.reactions` script, for which only the excitation energy corresponding to the current folder needs to be specified.

The trajectory directory structure is as follows:

```text
<energy>/
    <trajectory>/
        ├── geom.xyz
        ├── geom_$num$.xyz
        ├── movie.xyz
        ├── movie_$num$.xyz
        ├── reactions.dat
        ├── angular_momentum.dat
        ├── rotational_energy.dat
        ├── translational_energy.dat
        └── vibrational_energy.dat
```

For subsequent analyses, it is important to create an `analysis.pkl` file. This procedure helps filter out highly stretched and strongly off-equilibrium geometries by clustering structures that share the same SMILES descriptor. The clustering is performed using `sort.py`, with `r.sort` serving as the corresponding submission script. The excitation-energy folders to be included in the analysis must be specified in this script. In most cases, all available excitation energies should be selected.

Once the sorting procedure is complete, a `geoms` folder is created. This folder contains multiple files named according to their corresponding SMILES descriptors (i.e., `SMILES.dat`). Each file contains the structures assigned to that particular SMILES descriptor.

Next, a representative structure is determined for each cluster, with the medoid used as the representative. This step is performed using `medoids.py` and the `r.medoids` submission script. No additional user input is required, provided that the `geoms` folder containing the clustered structures is present. This procedure creates a `<medoids>` folder containing files named according to their corresponding SMILES descriptors. Unlike the files in `geoms`, each file in `<medoids>` contains only one structure: the medoid of the corresponding cluster.

The medoid structures are then optimized using the `opt.py` script and the `r.opt` submission script. No additional user input is required. This final step produces the `analysis.pkl` file, which maps the SMILES descriptors obtained from the trajectories to those assigned after geometry optimization. In addition, an `<opts>` folder is created. This folder contains files named according to their corresponding SMILES descriptors, with each file containing the optimized structure of the respective medoid. Overall, this procedure helps filter out highly distorted and off-equilibrium structures.

To extend existing trajectories by increasing the number of simulation steps, use the `-r` restart flag. When restarting existing trajectories, there is no need to specify the number of new trajectories using the `-n` flag. The script automatically extracts the final frame from `movie.xyz` and uses it to create `geom1.xyz`. Similarly, the final frame of `movie1.xyz` is used to create `geom2.xyz`, and the procedure continues analogously for subsequent restart files. No separate command is required because the script automatically detects the existing trajectory files and determines the appropriate restart index.

# Analysis

After all of the preceding steps have been completed, further analyses can be performed. These include analyzing the generated fragments, their energy partitioning and angular-momentum distributions, reaction pathways, species populations, and reaction rate constants.

The occurrence of four general channels—intact furan, isomerization, fragmentation, and fragmentation involving H/H₂ loss—as well as the number of fragments produced is visualized using `frag_dist.py`. Plotting parameters can be modified directly in the marked plotting sections of the code.

A more detailed analysis of the species and fragmentation channels is provided by `frag_dist_total.py`. Plotting parameters can likewise be modified in the marked sections of the code.

Energy partitioning and angular-momentum distributions are analyzed using `distributions.py`, with the energy values are specified through command-line arguments.

Reaction pathways are visualized as a node graph using `plot_reactions.py`. This reaction-pathway analysis tool represents chemical reaction networks using graph theory, where nodes represent chemical species and edges represent reactions between them. The analysis can be performed for one or more excitation-energy directories and provides several options for filtering and customizing the resulting network. Chemical species are represented by nodes, while chemical reactions are represented by edges. The number of reactions is displayed on each edge, and the edge width scales with the square root of the number of reactions. For species where assigning bond orders and bond connectivity was difficult, the node frames are highlighted in red. The SMILES descriptions of these species end with “NeighbourList", indicating that bond orders could not be reliably assigned and that all bonds are therefore represented with a bond order of one.

Users can specify the number of fragmentation events, exclude pathways involving hydrogen loss, and set the minimum number of occurrences required for a reaction to be displayed. The visualization can also be customized by adjusting the molecular image frame size, overall figure dimensions, and label offset. Molecular nodes can optionally be displayed using compact labels, with the corresponding label-to-SMILES legend printed to standard output.

For example, the following command visualizes all reaction pathways at an excitation energy of 8 eV and displays compact labels, with the corresponding legend printed to standard output:

```bash
py plot_reactions.py -e 8 -l
```

Use the `-f` flag to specify the minimum number of fragmentation events. To visualize only isomerization, use:

```bash
py plot_reactions.py -e 8 -f 0
```

To visualize reactions associated with fragmentation events while excluding reactions involving hydrogen loss (`-n`), and to display only reactions occurring at least four times (`-m`), use:

```bash
py plot_reactions.py -e 8 -f 1 -n -m 4
```

For more information about these flags, including options controlling graphical aspects of the visualization, use the `-h` flag.

Species populations and chemical rate constants in the corresponding microcanonical ensemble are calculated using `kinetics.py`. This script also creates a `<populations>` folder containing files with the populations of all species at each excitation energy, allowing the results to be visualized externally as needed.

Similarly, a `<rates>` folder is created when the `-r` flag is used. The files in this folder contain the rate constants for all identified chemical reactions, expressed in ps⁻¹.

For example, to generate the population data using an effective time interval of 10 fs between saved frames and 3000 simulation steps, use:

```bash
py kinetics.py -e 7 -t 10 -s 3000
```

The effective timestep is used to define the time interval between consecutive saved frames and therefore ensures that the simulation time is visualized correctly.

To visualize only the eight most abundant populations (`-m`) while also calculating and recording the reaction rate constants, use:

```bash
py kinetics.py -e 7 -t 10 -s 3000 -m 8 -r
```

While calculating the rate constants the code will also visalizes the generator matrix and stochimotry matrix in grahpical repsentation. Moreover, you will also be able to obsere the comparsin between the original populations and those caluclated using the rate constant for your check.

For more information, use the `-h` flag.

To compare populations previously obtained using `kinetics.py` and stored in the `<populations>` folder, use the `populations.py` script. For example, to compare populations at excitation energies of 8, 9, and 10 eV using an effective timestep of 10 fs, use:

```bash
py populations.py -e 8 9 10 -t 10
```

The number of steps is determined automatically from the population files.

The `-m` flag can also be used to visualize only a specified number of the most abundant populations. A complete guide to the available options can be displayed using the `-h` flag.

Finally, conservation of the total energy and total angular momentum is visualized using `thermal.py`. The left column displays the total energy relative to its value in the first frame, while the right column displays the magnitude of the total angular momentum.

For example, to visualize trajectories at excitation energies of 7 and 8 eV using an effective time interval of 10 fs between saved frames, use:

```bash
py thermal.py -e 7 8 -t 10
```

# Tutorial

This section demonstrates a typical workflow by presenting the commands in the order in which they are generally used while reiterating the most important guidelines described above.

First, a molecular geometry file for the molecule to be thermally excited is required. Consider furan with an input geometry stored in `furan.xyz`. The supplied geometry does not need to be optimized beforehand because geometry optimization is performed automatically by the script.

Suppose we want to simulate furan at excitation energies of 7 and 8 eV, with 100 trajectories at each energy. The corresponding initial conditions can be generated using:

```bash
py get_sets.py -i furan.xyz -e 7 8 -n 100
```

For linear molecules, add the `-l` flag.

At this stage, the script generates the initial conditions for the trajectories. Each excitation energy has its own folder, which contains trajectory subfolders with their respective initial conditions stored in `geom.xyz`.

The trajectories must then be propagated. This is performed using `omol.py`, which is deployed on HPC and controlled by the `r.omol` submission script. The user specifies the required parameters in the `CONFIGURATION` section:

```bash
########## CONFIGURATION ##########

STAGE="10"

GEOM="geom.xyz" # Input file
OUTPUT="movie.xyz" # Output file

TIMESTEP=0.1 # fs
STEPS=250000

CHARGE=0
SPIN=1
INTERVAL=100 # Recording frame interval (stride)

DATASET="..."

###################################
```

`STAGE` represents the trajectory folder. Thus, for trajectory folder `<10>`, this variable should be set to `10`.

`GEOM` specifies the file containing the initial conditions, while `OUTPUT` defines the name of the trajectory file containing the saved simulation frames. For the initial simulation, the input file is typically named `geom.xyz`, while the trajectory output is stored in `movie.xyz`.

`TIMESTEP` specifies the integration timestep in femtoseconds, and `STEPS` specifies the total number of integration steps. With a timestep of 0.1 fs and 250,000 steps, the total simulation time is 25 ps.

`CHARGE` and `SPIN` specify the charge and spin state of the system.

`INTERVAL`, also referred to as the stride, determines how frequently simulation frames are recorded. In this example, the timestep is 0.1 fs and the stride is 100, so a trajectory frame is recorded every 10 fs of simulation time.

`DATASET` specifies the path to the dataset used by the eSEN machine-learned interatomic potential.

Once these parameters have been configured, the simulation is ready to be deployed on the computational cluster.

After the simulation is complete, the conservation of total energy and total angular momentum can be inspected using:

```bash
py thermal.py -e 7 8 -t 10
```

For the subsequent analysis, we need to identify all chemical species and calculate additional observables, such as energy partitioning and the angular momenta of the furan isomers and their fragments. This analysis is performed using the `reactions.py` script, which is deployed on the computational cluster using the `r.reactions` submission file. In the `CONFIGURATION` section of `r.reactions`, set the `ENERGY` variable to the excitation-energy folder that you want to analyze. Once the calculations are complete, the script automatically transfers the resulting files to the corresponding `<trajectory>` folders.

Moreover, the filtering procedure described in the second section must first be performed. Structures sharing the same SMILES descriptor are grouped into common files within the `<geoms>` folder using `sort.py`. The script is deployed on the computational cluster using the `r.sort` submission file, in which the excitation-energy folders to be processed must be specified. In this case, the excitation energies from 7 to 8 eV are selected using:

```bash
ENERGY==$(seq 8 -1 7)
```
Next, a medoid representing each cluster is calculated using `medoids.py`, which is submitted using `r.medoids`. The resulting medoid structures are stored in the `<medoids>` folder.

The medoids are subsequently optimized using `opt.py` and `r.opt`. This procedure generates the `analysis.pkl` file, which maps the SMILES descriptors of the original trajectory clusters to the SMILES descriptors obtained after optimization. This mapping helps filter highly distorted and off-equilibrium geometries. The `analysis.pkl` file also contains frequencies obtained from the frequency analysis, which can be useful for identifying and post-processing structures that did not optimize correctly.

Next, fragmentation events and the number of fragments produced can be analyzed using:

```bash
py frag_dist.py
```

A more detailed analysis of the generated species, isomers, and fragmentation channels can be performed using:

```bash
py frag_dist_total.py
```

Next, the energy-partitioning and angular-momentum distributions of the fragments can be visualized for the excitation-energy folders of interest. For example, for 7 and 8 eV:

```bash
py distributions.py 7 8
```

Reaction pathways are analyzed using the `plot_reactions.py` script. To visualize the reaction pathways at 7 eV, use:

```bash
py plot_reactions.py -e 7
```

Reaction pathways at higher excitation energies may become more complex, making it useful to filter the displayed reactions. For example, to display only pathways before fragmentation at 8 eV, use:

```bash
py plot_reactions.py -e 8 -f 0
```

To display pathways involving at least one fragmentation event while excluding reactions associated with hydrogen loss, use:

```bash
py plot_reactions.py -e 8 -f 1 -n
```

You may also want to display only the most frequently occurring reaction pathways. For example, to display reactions with more than five occurrences, use:

```bash
py plot_reactions.py -e 8 -r 5
```

For easier interpretation of complex networks, compact labels can be enabled using `-l`. The label offset and molecular image frame size can also be customized using the corresponding graphical options. For example:

```bash
py plot_reactions.py -e 8 -l -o 20 -i 240
```

For reaction analysis, it is also useful to visualize species populations. First, population data must be generated and stored in the `<populations>` folder. For example, to compute the populations at an excitation energy of 7 eV using an effective timestep of 10 fs over a total simulation steps of 2500, use:

```bash
py kinetics.py -e 7 -t 10 -s 2500
```

This script also visualizes the populations at the selected excitation energy, providing a convenient check of the generated data.

Reaction rate constants can additionally be calculated and recorded using the `-r` flag:

```bash
py kinetics.py -e 7 -t 10 -s 2500 -r
```

The rate constants for the identified reactions are then stored in the `<rates>` folder within the corresponding excitation-energy directory.

Once the `<populations>` folders have been generated, populations at different excitation energies can be compared. For example, populations at 7 and 8 eV can be compared using:

```bash
py populations.py -e 7 8 -t 10
```

For systems containing a large number of species, only a specified number of the most populated species can be visualized:

```bash
py populations.py -e 7 8 -t 10 -m 8
```

In this example, the eight most abundant chemical species are displayed for each excitation energy.

After analyzing the results, it may become necessary either to add more trajectories or to extend the existing trajectories.

To extend the existing trajectories in, for example, the 8 eV folder, add the `-r` flag to `get_sets.py`:

```bash
py get_sets.py -i furan.xyz -e 8 -r
```

This generates a new initial-condition file, `geom1.xyz`, for each existing trajectory using the final frame of `movie.xyz`. The procedure can be repeated for subsequent restarts. For example, if the most recent trajectory file is `movie9.xyz`, the script automatically generates `geom10.xyz`.

To add additional trajectories, use:

```bash
py get_sets.py -i furan.xyz -e 8 -n 200
```

When `-n` is used, the script creates the specified number of additional trajectories—in this example, 200.

When propagating restarted trajectories, remember to update the corresponding settings in `r.omol`. In particular, `GEOM` and `OUTPUT` must refer to the appropriate restart files. For example, when restarting from `movie.xyz` using the newly generated `geom1.xyz`, set:

```bash
GEOM="geom1.xyz" # Input file
OUTPUT="movie2.xyz" # Output file
```

and proceed analogously for subsequent restarts. All previously used analysis scripts will automatically detect and load all available `movie*.xyz` files. However, the analyses must be rerun to incorporate the newly generated trajectory data.

This represents the typical workflow for carrying out and analyzing these simulations. All scripts controlled through command-line flags also provide a `-h` flag that displays the available arguments and additional usage information.

