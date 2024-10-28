import matplotlib.pyplot as plt
from matplotlib import cm, colors
from Bio import Phylo
import numpy as np
from msa_utils import *

# plot output of phylop
def plot_phylop_scores(phylop_data, output_plot_file):
    """
    Plots PhyloP scores vs. position and saves the plot to a file.

    Parameters:
    - phylop_data: List of dictionaries containing PhyloP data, as output by read_phylop_output_fixed_step().
    - output_plot_file: Path to save the plot as an image file.
    """
    # Extract positions and scores from the data
    positions = [entry["position"] for entry in phylop_data]
    scores = [entry["scale"] for entry in phylop_data]

    # Create the plot
    plt.figure(figsize=(10, 6))
    plt.plot(positions, scores, label="PhyloP Score", color="b", linewidth=0.8)

    # Add labels and title
    plt.xlabel("Position")
    plt.ylabel("PhyloP Score")
    plt.title("PhyloP Scores vs. Position")

    # Add a grid for better readability
    plt.grid(True)

    # Save the plot to a file
    plt.savefig(output_plot_file, dpi=300)

    # Close the plot to avoid displaying it in environments that support inline plotting
    plt.close()


def get_tree_y_positions(tree):
    """Return a dictionary of leaf names and their y-positions in the tree plot."""
    leaves = tree.get_terminals()
    y_positions = {leaf.name: idx for idx, leaf in enumerate(leaves)}
    return y_positions


def color_tree(tree, subtree_name=None, values=None, msa=None,
               cmap_name='coolwarm', output_file="out_tree.png"):
    """
    Generalized function to color a tree.
    - Either provide a subtree name to color that subtree.
    - Or provide a vector of numerical values to color branches with a heatmap.

    Parameters:
    - tree: A Bio.Phylo tree object.
    - output_file: Path to save the output image.
    - subtree_name: Name of the internal node whose subtree will be colored.
    - values: A list or array of numerical values (one for each branch).
    - msa: A Bio.Align.MultipleSeqAlignment object containing the alignment (optional).
           If provided, it will be displayed with rows aligned to tree leaf nodes.
    - cmap_name: Name of the colormap to use for heatmap coloring.
    - output_file: Name of file to save the figure.
    """
    # Prepare color map if using numerical values
    cmap = plt.get_cmap(cmap_name)

    # Case 1: Color by subtree
    if subtree_name:
        # Function to color branches of a clade and its subclades
        def color_clades(clade, color):
            for subclade in clade.clades:
                subclade.color = color  # Add color to each clade
                color_clades(subclade, color)

        # Find the internal node by name
        found_clade = None
        for clade in tree.find_clades():
            if clade.name == subtree_name:
                found_clade = clade
                break

        if found_clade is None:
            raise ValueError(f"Internal node {subtree_name} not found in the tree.")

        # Color the entire tree (default color black)
        for clade in tree.get_terminals() + tree.get_nonterminals():
            clade.color = 'black'  # Set the default color to black for all nodes

        # Color the subtree of the found internal node in red
        color_clades(found_clade, 'red')

    # Case 2: Color by numerical values (like a heatmap)
    elif values is not None:
        values_list = list(values.values())

        # Normalize and get the colormap
        norm = colors.Normalize(vmin=np.min(values_list), vmax=np.max(values_list))
        cmap = plt.get_cmap(cmap_name)

        values_dict_simple = {k.name: values[k] for k in values.keys()}
        # Map each clade to its color based on the normalized values
        for clade in tree.get_nonterminals() + tree.get_terminals():
            if clade.name in values_dict_simple:
                color_value = norm(values_dict_simple[clade.name])
                branch_color = cmap(color_value)
                # Assign color to clade
                clade.color = colors.rgb2hex(branch_color)


        # Normalize the values to [0,1] for the color map
#        norm = colors.Normalize(vmin=np.min(values), vmax=np.max(values))
#        branch_colors = cmap(norm(values))

        # Assign colors to each branch based on the values
#        for clade, branch_color in zip(tree.get_nonterminals() + tree.get_terminals(), branch_colors):
#            # Convert the RGB tuple (from colormap) to a hex string for matplotlib
#            clade.color = colors.rgb2hex(branch_color)

    else:
        raise ValueError("Either subtree_name or values must be provided.")

        # Custom label function to return the name for terminal nodes, and set leaf color

    def leaf_labels(clade):
        if clade.is_terminal():
            # Return the clade's name if terminal (leaf) and adjust color accordingly
            return clade.name
        return None

    # Draw the tree
    fig, ax_tree = plt.subplots(figsize=(10, 10))
#    if msa is not None:
    print("Adding MSA!!!")
    ax_msa = fig.add_axes([0.6, 0.1, 0.3, 0.8])  # Add MSA subplot

    Phylo.draw(tree, label_func=leaf_labels, do_show=False, axes=ax_tree)

#    fig = plt.figure(figsize=(10, 10))
#    axes = fig.add_subplot(1, 1, 1)
    # Draw the tree with colored branches
#    Phylo.draw(tree, axes=axes, label_func=leaf_labels,  do_show=False)

    # If coloring by numerical values, add a color bar
    if values is not None:
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])  # Set a dummy array for the colorbar
        cbar = plt.colorbar(sm, ax=ax_tree)
        cbar.set_label('Branch Values')

    for text in ax_tree.texts:  # axes.texts contains all the text objects
        leaf_name = text.get_text()
        for clade in tree.get_terminals():
            if clade.name == leaf_name.strip():
                # Set the color of the text based on the clade's assigned color
                if isinstance(clade.color, str):  # color names (red, green ...)
                    text.set_color(clade.color)
                else:  # color numbers 0-255
                    text.set_color((clade.color.red / 255, clade.color.green / 255, clade.color.blue / 255))

    print("BEFORE IF MSA")
    if msa is not None:  # Display also msa sequences
        print("Adding MSA!!!")
#        ax_msa = fig.add_axes([0.6, 0.1, 0.3, 0.8])  # Add MSA subplot

        # Get y-positions of leaves to align MSA correctly
        y_positions = get_tree_y_positions(tree)
        print("y_positions before: ", y_positions)
        y_positions = {key.rstrip('0123456789'):y_positions[key] for key in y_positions}
        print("y_positions after: ", y_positions)

        # Define nucleotide colors
        nucleotide_colors = {"A": "green", "C": "blue", "G": "orange", "T": "red", "-": "lightgrey"}

        # Collect leaf names in the order they appear in the tree drawing
        leaves = tree.get_terminals()
        leaf_names = [leaf.name for leaf in leaves]

        # Extract species names from the tree
        tree_species = {leaf.name for leaf in tree.get_terminals()}

        # Filter MSA blocks by tree species
        print("Tree Species: ", tree_species)
        print("Length msa before:")
###        msa = filter_msa_blocks_by_species(msa, tree_species)

        # Plot the MSA blocks, offsetting each block sequentially
        x_offset = 0
        for block in msa:
#            print("Filtered MSA Block: ")
#            print(block)
            for record in block:
                species_name = record.id.split(".")[0].rstrip('0123456789')  # remove genome version
#                print("Trying to find ", species_name)
                if species_name in y_positions:
                    print("Add patch!!! ", str(record.seq))
                    print("species_name: ", species_name, " record.id: ", record.id, " y.positions: ", y_positions)
                    y = y_positions[species_name]
                    for x, nucleotide in enumerate(str(record.seq)):
                        color = nucleotide_colors.get(nucleotide, "black")
                        ax_msa.add_patch(plt.Rectangle((x + x_offset, y), 1, 1, color=color))
            x_offset += block.get_alignment_length()
            break  # Take one block!! Temporary!!

        # Formatting
#        print("Got y positions: ", y_positions)
#        print("Got x_offset: ", x_offset)
        ax_msa.set_xlim(0, x_offset)
        ax_msa.set_ylim(0, len(y_positions))
        ax_msa.axis("off")  # Turn off axis for the MSA
        ax_msa.set_title("MSA")

#        if isinstance(msa, list): # multiple blocks !!! take first !!
#            use_msa = msa[0]
#        else:
#            use_msa = msa
#        # Create a dictionary mapping leaf names to their MSA sequences (MSA can come in blocks!!!)
#        msa_dict = {record.id: str(record.seq) for record in use_msa}

#        # Create an alignment matrix to display the sequences alongside the tree
#        msa_matrix = [list(msa_dict[leaf_name]) if leaf_name in msa_dict else [''] * use_msa.get_alignment_length() for
#                      leaf_name in leaf_names]

        # Convert MSA matrix to numpy for easier display
 #       msa_array = np.array(msa_matrix)

        # Create a new axes for the MSA display
#        ax_msa = fig.add_axes([0.8, 0.1, 0.15, 0.8])  # Adjust position and size as needed
#        ax_msa.set_title("MSA")
#        # Display the MSA matrix as an image
#        # Color mapping for nucleotides or amino acids
#        msa_cmap = plt.cm.Pastel1
#        ax_msa.imshow(msa_array == 'A', cmap=msa_cmap, aspect='auto',
#                      interpolation='none')  # Example with 'A'; update to generalize
#        # Adjust the axis
#        ax_msa.set_xticks(range(use_msa.get_alignment_length()))
#        ax_msa.set_yticks(range(len(leaf_names)))
#        ax_msa.set_yticklabels(leaf_names, fontsize=8)
#        ax_msa.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)

        # NEW: Create MSA heatmap with text overlay
 #       ax_msa = fig.add_axes([0.75, 0.1, 0.15, 0.8])
 #       msa_array = np.array([list(rec.seq) for rec in use_msa])
 #       msa_data = np.vectorize(lambda x: ord(x))(msa_array)  # Convert characters to int for color mapping

        # Display MSA heatmap
 #       im = ax_msa.imshow(msa_data, aspect='auto', cmap='coolwarm')
 #       plt.colorbar(im, ax=ax_msa, orientation="vertical", label="MSA")

        # Overlay text of each character in MSA
#        for i in range(msa_data.shape[0]):
#            for j in range(msa_data.shape[1]):
#                ax_msa.text(j, i, msa_array[i, j], ha='center', va='center', color='black', fontsize=6)

    # Save the plot to a file
    print("SAVING TREE+MSA!!!")
    plt.savefig(output_file)
    plt.close()


def subtrees_to_color_vector(tree, subtrees):
    """
    Creates a color vector for branches in `tree` based on a list of `subtrees`.

    Parameters:
    - tree: A Bio.Phylo tree object containing all species.
    - subtrees: A list of Bio.Phylo tree objects, each representing a non-overlapping subtree of `tree`.

    Returns:
    - A list where each integer represents a color code for a clade:
      0 for branches not in any subtree, 1 for the first subtree, 2 for the second, etc.
    """
    # Initialize dictionary to hold color codes for each clade in the tree
    values_dict = {}

    # Assign unique color to each subtree in the input list
    for color_index, subtree in enumerate(subtrees, start=1):
        subtree_species = {leaf.name for leaf in subtree.get_terminals()}
        subtree_nodes = {clade.name for clade in subtree.find_clades() if clade.name}

        for clade in tree.find_clades():
#            clade_leafs = set(clade.name.split("-")).issubset(subtree_species)
            # Assign the color only if the clade is a terminal node and in the current subtree
#            if clade.name in subtree_nodes:  # clade.is_terminal() and
            if set(clade.name.split("-")).issubset(subtree_species):
                values_dict[clade] = color_index

    # Assign color 0 to all branches not part of any subtree
    for clade in tree.find_clades():
        if clade not in values_dict:
            values_dict[clade] = 0

    # Convert values_dict to a list vector for the color_tree function
#    values_vector = [values_dict[clade] for clade in tree.find_clades()]

    return values_dict  # values_vector

