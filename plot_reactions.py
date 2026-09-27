"""
Analyze molecular reaction trajectories using graph-based reaction pathways.

The script loads species assignments from trajectory-specific ``reactions.dat``
files together with a species mapping stored in ``analysis.pkl``. It smooths
frame-by-frame species trajectories, compresses repeated states into reaction
sequences, filters selected pathways, and builds a directed reaction network.
Repeated transformations are counted and represented as weighted graph edges.

RDKit is used to generate two-dimensional molecular depictions for graph nodes,
while NetworkX and Graphviz provide the reaction-network representation and
layout. The final reaction pathway is written to ``Reaction_pathways.png`` and
summary information, including optional species aliases and node populations,
is printed to standard output.

Command-line arguments control the analyzed energy folders, fragmentation
filtering, treatment of hydrogen-loss pathways, graphical layout parameters,
node labels, molecular image size, and the minimum number of reactions required
for an edge to be displayed.
"""

	# %%
# Numerical operations, plotting, filesystem handling, graph construction,
# molecular depictions, serialization, tabular output, and CLI parsing.
import numpy as np
import matplotlib.pyplot as plt
import os
import networkx as nx
from collections import Counter
import networkx as nx
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import Draw
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
from rdkit.Chem.Draw import IPythonConsole
# Enable RDKit 3D rendering support when this script is executed interactively.
IPythonConsole.ipython_3d = True
from networkx.drawing.nx_agraph import graphviz_layout
from matplotlib.offsetbox import OffsetImage, AnnotationBbox, TextArea
from PIL import ImageDraw, ImageFont
import pickle
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from PIL import Image
from io import BytesIO
import pandas as pd
import argparse




# Configure command-line options for pathway filtering and visualization.
parser = argparse.ArgumentParser(
                prog="Node graph theory-based analysis of chemical reaction pathways.",
                description="""This reaction pathway analysis tool visualizes chemical reaction pathways using graph theory. Each node represents a chemical species,
                while the edges between nodes represent chemical reactions. Each edge displays the number of reactions, with its width scaled proportionally to the 
                square root of the number of reactions. The analysis provides several options for custom filtering. You can specify the minimum number of fragmentation events,
                choose whether to display reaction pathways involving hydrogen loss, and set the minimum number of reactions between two chemical species required for an edge to
                be displayed. Furthermore, several graphical parameters can be adjusted to customize the resulting visualization."""
            )



# One or more energy-directory names must be supplied.
parser.add_argument("-e", "--energies", nargs="+", required=True, type=str, help="Enter one or more energy directories to be used for the analysis.")

# Optional threshold controlling how strongly fragmented pathways are retained.
parser.add_argument("-f", "--fragmenation", type=int, help="Enter the minimum number of fragmentation events required for a pathway to be displayed.")

# Toggle the pathway-selection branch associated with terminal hydrogen species.
parser.add_argument("-n", "--no_hydrogen_loss", action="store_true", help="Use this flag if you do not want to display reaction pathways involving hydrogen loss.")

# Vertical offset used when positioning textual node aliases.
parser.add_argument("-o", "--offset", nargs="?", type=int, default=15, help="Enter the offset between the label and the image frame. The default is 15.")

# Pixel dimensions of each RDKit molecular drawing.
parser.add_argument("-s", "--frame_size", nargs="?", type=int, default=280, help="Enter the size of the molecular image frame in pixels. The default is 280.")

# Request display of compact S-number labels beneath molecular nodes.
parser.add_argument("-l", "--display_label", action="store_true", help="Use this flag to display labels for molecular nodes. A legend linking the labels to their corresponding SMILES will be printed to standard output.")

# Overall Matplotlib figure width and height.
parser.add_argument("-i", "--image_size", nargs=2, default=[12, 8], type=int, help="Enter the width and height of the plot figure. The default is 12 8.")

# Minimum repeated reaction count required before an edge is plotted.
parser.add_argument("-r", "--min_reactions", nargs="?", default=0, type=int, help="Enter the minimum number of reactions required between two chemical species for the corresponding edge to be displayed. The default is 1.")





# Parse the command line once all supported options have been registered.
args = parser.parse_args()


# Store parsed command-line arguments in local variables used below.
set_folders = args.energies #sorted([folder for folder in os.listdir() if folder.isdigit()], key=int)
min_frag_events = args.fragmenation
hydrogen_loss = args.no_hydrogen_loss
offset = args.offset
frame_size = args.frame_size
display_label = args.display_label
image_size = args.image_size
min_reactions = args.min_reactions


# %%
# Collect frame-resolved species assignments for each energy and trajectory.
reactions = {}


# Load the mapping used to convert stored molecular identifiers into species names.
with open("analysis.pkl", "rb") as f:
    data = pickle.load(f)
    
    # Normalize the special NeighborList hydrogen identifier to the alias mapping used below.
    data["[HH]NeighborList"] = ["[H]"]


# Process each requested energy set independently.
for set_ in set_folders :
    reactions[set_] = {}

    # Each entry inside an energy directory is treated as a candidate trajectory.
    trajectories = os.listdir(set_)

    for trajectory in trajectories :

        # Ignore trajectory directories that do not contain the expected reaction file.
        if os.path.isfile(f"{set_}/{trajectory}/reactions.dat") :
            reactions[set_][trajectory] = []

            with open(f"{set_}/{trajectory}/reactions.dat") as f :
                for line in f.readlines() :

                    # Skip the first tab-separated field (frame index) and retain the listed species.
                    elements = list(map(str, line.rstrip("\n").split("\t")[1:]))

                    # Translate each stored molecular identifier through the preloaded mapping.
                    reactions[set_][trajectory].append([data[mol][0] for mol in elements])

            print(f"File {set_}/{trajectory}/reactions.dat was loaded.")



# Smooths the species trajectory
# Replace inconsistent transient assignments with the preceding accepted state.
names = {}

# Collect unique canonical species names before assigning short display aliases.
all_names = set([i[0] for i in data.values()])
for i, line in enumerate(sorted(all_names)) :
    # Alias ordering is deterministic because the names are sorted first.
    names[line] = "S" + str(i)



# Number of frames considered by the local consistency check below.
steps_ahead = 1
trajs = {}
for energy, set_ in reactions.items() :

    trajs[energy] = {}
    for trajectory, reaction in set_.items() :
        # Work on the frame sequence associated with the current trajectory.
        frames = [a for a in reaction]

        # Seed the smoothed trajectory with the first observed state.
        traj = [frames[0]]
        for i in range(1, len(frames)) :
            if i < len(frames) :
                
                # Keep the current frame when its species set agrees with the selected look-ahead window.
                if set(list(set(np.concatenate(frames[i:i+steps_ahead])))) == set(np.array(frames[i])) :
                    traj.append(frames[i])
                else :
                    # Otherwise carry forward the previously accepted state.
                    traj.append(traj[i-1])
                    
                # Every specie with "." which means that RDkit defines that as framgnets
                # A dot in an RDKit SMILES denotes disconnected components; suppress such assignments here.
                if any("." in element for element in traj[i]) : #and len(traj[i-1]) == len(traj[i]) :
                    traj[i] = traj[i-1]                    

        # Store the fully smoothed frame sequence for later compression.
        trajs[energy][trajectory] = traj


# Shortens the blocks of same species
# Compress consecutive identical frames into discrete state-change events.
sequences = {}
for energy, set_ in trajs.items() :

    sequences[energy] = {}
    for trajectory, traj in set_.items() :

        frames = traj

        # Record the initial state together with the frame at which it begins.
        products = [(0, frames[0])]

        for i in range(1, len(frames)) :
            # A new product/state is recorded only when the frame composition changes.
            if frames[i-1] != frames[i] :
                products.append((i, frames[i]))
        
        to_remove = []
        for j, product in enumerate(products[1:-1], 1) :
            # Inspect a three-state window to identify short-lived fragmentation excursions.
            prev, current, next = products[j-1][1], products[j][1], products[j+1][1]

            # Filtres out the reactions such A + B -> C + D + E -> A + B
            if len(current) > 1 :
                
                # Mark an expanded intermediate for removal when the surrounding states have equal size.
                if len(prev) == len(next) and len(prev) < len(current) and len(next) < len(current) :
                    to_remove.append(j)

        # Delete from the end so earlier list indices remain valid.
        for k in sorted(to_remove, reverse=True):
            products.pop(k)
        


        # Save the compressed state-change sequence for this trajectory.
        sequences[energy][trajectory] = products



# Removes the blocs between same species
# Remove loops through previously visited states and apply pathway-selection filters.
all_structures = []

reaction_path = {}
for energy, set_ in sequences.items() :

    reaction_path[energy] = {}
    for trajectory, sequence in set_.items() :

        result = []
        
        # Count species present in the terminal state for the final population summary.
        all_structures.extend([names[i] for i in sequence[-1][1]])
        
        # Standard branch: retain trajectories without requiring a terminal hydrogen species.
        if not hydrogen_loss :
             
            for item in sequence :
                # Convert the species list into a tuple so repeated states can be compared directly.
                key = tuple(item[1])

                # Build the sequence of states already accepted into the current reaction path.
                keys = [tuple(x[1]) for x in result]

                # Revisiting a previous state indicates a loop in the trajectory.
                if key in keys :
                    idx = keys.index(key)
                    # Truncate the path back to the first occurrence of the revisited state.
                    result = result[:idx + 1]
                else :
                    result.append(item)
                    
        else :
         
            # In the alternate branch, only terminal states containing H2 or H are processed.
            if "[H][H]" in sequence[-1][1] or "[H]" in sequence[-1][1] : 
            
                for item in sequence :
                    key = tuple(item[1])

                    keys = [tuple(x[1]) for x in result]

                    if key in keys :
                        idx = keys.index(key)
                        result = result[:idx + 1]
                    else:
                        result.append(item)
         
            
      
        
        # Without a fragmentation threshold, retain the loop-reduced sequence as-is.
        if min_frag_events == None :
            reaction_path[energy][trajectory] = result 
            
        else :
            
            # A zero threshold keeps only states consisting of a single species.
            if min_frag_events == 0 :
                reaction_path[energy][trajectory] = [i for i in result if len(i[1]) == 1]

            else :
                
                reaction_path[energy][trajectory] = []
        
                
                # Once the requested fragmentation level is crossed, subsequent states are retained.
                breakpoint = False
                for i in range(len(result)) :
                    step = result[i]

                    # Detect the first state whose fragment count exceeds the requested threshold.
                    if len(step[1]) > min_frag_events :
                        # Include the immediately preceding state to preserve the transition into fragmentation.
                        reaction_path[energy][trajectory].append(result[i-1])
                        breakpoint = True
                        
                        
                    # After the threshold has been reached, append the remaining pathway states.
                    if breakpoint :         
                        reaction_path[energy][trajectory].append(step)
                    
                    
            
            
print("-"*60)

for energy, traj in reaction_path.items() :

    for n, elements in traj.items() :
        seq = []
        for element in elements :
            # Flatten the selected state contents; this preserves the original reporting structure.
            seq.extend(element[1])


# Build a multigraph that preserves every observed fragment-to-fragment reaction event.
G = nx.MultiDiGraph()
for energy, set_ in reaction_path.items() :
    for trajectory, sequence in set_.items() :

        # Compare each pair of consecutive reaction states.
        for (t1, s1), (t2, s2) in zip(sequence[:-1], sequence[1:]):
            # Count multiplicities so repeated identical fragments are handled correctly.
            c1 = Counter(s1)
            c2 = Counter(s2)

            # Species present before but not after the transition are reactant-side losses.
            lost = list((c1 - c2).elements())
            # Species appearing after the transition are product-side gains.
            gained = list((c2 - c1).elements())

            #G.add_edge("A", "B")
            
            #if len(lost) <= len(gained) :
            
            # Connect every disappearing fragment to every newly appearing fragment.
            for frag_1 in lost :
                #if "NeighborList" not in frag_1 : 
                G.add_node(frag_1)

                for frag_2 in gained :
                    #if "NeighborList" not in frag_2 :
                    G.add_node(frag_2)
                        
                    #if "NeighborList" not in frag_1 and "NeighborList" not in frag_2 :
                    # Preserve the originating energy/trajectory on each individual multigraph edge.
                    G.add_edge(frag_1, frag_2, label=f"{energy}-{trajectory}")
                    
         
                   



# Collapse repeated multigraph edges into a weighted directed graph for plotting.
H = nx.DiGraph()

# Creates edges

for i, node_1 in enumerate(G.nodes()) :
    for j, node_2 in enumerate(G.nodes()) :

        # Count how many individual observations produced this directed transformation.
        edge_number = G.number_of_edges(node_1, node_2)
        # Exclude self-pairs and node pairs with no observed directed reaction.
        if i != j and edge_number != 0 :
            # Retain only transformations occurring more often than the configured threshold.
            if edge_number > min_reactions :  
            	# Store the reaction count for both line-width scaling and the displayed edge label.
            	H.add_edge(node_1, node_2, width=edge_number, label=edge_number)


# Compute a hierarchical layout suitable for directed reaction pathways.
pos = graphviz_layout(H, prog="dot")



# Draw
# Create the plotting canvas and render each molecular species as an RDKit image.
fg, ax = plt.subplots(figsize=(image_size[0], image_size[1]))


# Generate and place the molecular depiction associated with each graph node.
for node in H.nodes() : 
    
    #print(node)
    # NeighborList-derived species are handled separately and outlined in red.
    if "NeighborList" in node :
        # Remove the marker appended during fallback bond detection before RDKit parsing.
        node_smile = node.rstrip("NeighborList")
       
        # Parse without sanitization because fallback structures may not satisfy normal valence rules.
        mol = Chem.MolFromSmiles(node_smile, sanitize=False)
        
        color = "red"
        
        for atom in mol.GetAtoms():
            # Prevent RDKit from automatically adding implicit hydrogens to fallback structures.
            atom.SetNoImplicit(True)

    else :
        params = Chem.SmilesParserParams()
        # Preserve explicit hydrogen atoms when parsing standard species.
        params.removeHs = False
        # Avoid sanitization so the visualization follows the stored SMILES representation directly.
        params.sanitize = False
       
        mol = Chem.MolFromSmiles(node, params)
       
        color = "black"
        
    
    # Use RDKit's 2D coordinate generator
    # Atomic/molecular hydrogen receives special coordinate and line-width handling.
    if node == "[H][H]" or node == "[H]" : 
        rdDepictor.SetPreferCoordGen(False)
    else :
        rdDepictor.SetPreferCoordGen(True)
        
    # Generate two-dimensional coordinates before rasterizing the molecule.
    rdDepictor.Compute2DCoords(mol)

    # Render each molecule to an in-memory Cairo PNG of the requested size.
    drawer = rdMolDraw2D.MolDraw2DCairo(frame_size, frame_size) # 280
    opts = drawer.drawOptions()
    
    if node == "[H][H]" or node == "[H]" : 
        opts.bondLineWidth = 2.5
    else :
        opts.bondLineWidth = 1.5

    rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol)
    drawer.FinishDrawing()

    # Retrieve the rendered PNG bytes and convert them into a PIL image.
    png = drawer.GetDrawingText()
    img = Image.open(BytesIO(png))

    # Scale the molecular bitmap before embedding it at the graph-node position.
    imagebox = OffsetImage(img, zoom=0.25)

  
    # Place the molecular image inside a framed annotation centered on the node coordinate.
    ab = AnnotationBbox(imagebox, pos[node], frameon=True)
    
    # Border color encodes fallback structures (red), terminal nodes (green), or standard nodes (black).
    ab.patch.set_edgecolor(color)
    ax.add_artist(ab)



# Scale edge widths by the square root of the observed reaction count.
widths = [np.sqrt(H[u][v]["width"]) for u, v in H.edges()]

# Draw directed arrows between molecular nodes using the frequency-scaled widths.
edges = nx.draw_networkx_edges(
    H,
    pos,
    ax=ax,
    arrows=True,
    arrowstyle="-|>",
    arrowsize=10,
    min_source_margin=30,
    min_target_margin=30,
    alpha=0.4,
    width=widths
    #connectionstyle="arc3,rad=0.05"
)

# Raise reaction arrows above lower-priority plot elements.
for e in edges :
    e.set_zorder(10)

# Extract reaction counts for display on the corresponding directed edges.
edge_labels = nx.get_edge_attributes(H, "label")



# Draw frequency labels slightly before the midpoint of each edge.
edge_l = nx.draw_networkx_edge_labels(
    H,
    pos,
    edge_labels=edge_labels,
    label_pos=0.4,
    bbox=dict(
        facecolor="none",
        edgecolor="none")
    )

# Make labels appear on top
for text in edge_l.values():
    text.set_zorder(9)


# Optionally display compact species aliases beneath their molecular structures.
if display_label :
    # Move labels downward
    # Shift label coordinates downward by the user-configurable offset.
    label_pos = {
        node : (x, y - offset)
        for node, (x, y) in pos.items()
    }

    node_labels = nx.draw_networkx_labels(
        H,
        label_pos,
        # Display the compact S-number alias rather than the full SMILES string.
        labels={node : names[node] for node in H.nodes()},
        font_size=14,
        font_weight="bold"
    )

    for text in node_labels.values():
        text.set_zorder(11)


    print("SMILES and their aliases:")
    # Report populations only for species that are present in the plotted network.
    for node in H.nodes() :
        print(names[node], node)


    print("-"*60)

# Count how often each final structure appears across the processed trajectories.
# Aggregate the previously collected terminal-species aliases.
strucures_data = Counter(all_structures)


final_structures = {}


for node in H.nodes() :
    node_name = names[node]
    
    # Use the observed terminal count when the node appears among final structures.
    if strucures_data.get(node_name) is not None :
        final_structures[node_name] = strucures_data[node_name]
    else :
        # Nodes never observed as terminal products receive a zero population.
        final_structures[node_name] = 0


# Convert final node populations into a table for terminal output.
df = pd.DataFrame({
    "Name": final_structures.keys(),
    "Value": final_structures.values()
})

# Number of side-by-side Name/Value groups used when formatting the output table.
n_groups = 1
# Determine the number of rows needed per table group.
n_rows = int(np.ceil(len(df) / n_groups))

parts = []
for i in range(n_groups):

    # Slice each group and reset its index before horizontal concatenation.
    part = df.iloc[i*n_rows:(i+1)*n_rows].reset_index(drop=True)
    part.columns = ["Name", "Value"]
    parts.append(part)

# Combine the groups horizontally into the final printable population table.
wide_df = pd.concat(parts, axis=1)

print("Population of nodes:")
# ``to_string`` prevents pandas from truncating rows or columns in terminal output.
print(wide_df.to_string())



# Finalize, save, and display the reaction-pathway figure.
# Hide numeric axes because graph coordinates have no direct physical meaning.
ax.set_axis_off()
# Reduce unnecessary whitespace before saving the final figure.
plt.tight_layout()
# Save the complete reaction network with a tightly cropped bounding box.
fg.savefig("Reaction_pathways.png", bbox_inches="tight")
plt.show()






