import matplotlib.pyplot as plt
from matplotlib import cm, colors
from Bio import Phylo
import numpy as np
from msa_utils import *
from matplotlib import gridspec
import matplotlib.patches as patches


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
               cmap_name='rainbow', output_file="out_tree.png"):  # coolwarm
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

    if msa is None:
        # Draw the tree
        fig, ax_tree = plt.subplots(figsize=(10, 10))
    else: # add msa axes
        # Set up the grid layout to place MSA on the right
        fig = plt.figure(figsize=(10, 8))
        gs = gridspec.GridSpec(1, 2, width_ratios=[4, 3])  # Adjust width ratios for desired spacing

        # Tree and colorbar axis on the left
        ax_tree = fig.add_subplot(gs[0])

        # MSA plot axis on the right
        ax_msa = fig.add_subplot(gs[1])

    Phylo.draw(tree, label_func=leaf_labels, do_show=False, axes=ax_tree)
    ax_tree.set_ylabel('')

    # If coloring by numerical values, add a color bar
    if values is not None:
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])  # Set a dummy array for the colorbar
    #    if msa is None:
        cbar = plt.colorbar(sm, ax=ax_tree, location='left')
    #    else:
    #        cbar = plt.colorbar(sm, ax=ax_colorbar)
        cbar.set_label('Rate')

    for text in ax_tree.texts:  # axes.texts contains all the text objects
        leaf_name = text.get_text()
        for clade in tree.get_terminals():
            if clade.name == leaf_name.strip():
                # Set the color of the text based on the clade's assigned color
                if isinstance(clade.color, str):  # color names (red, green ...)
                    text.set_color(clade.color)
                else:  # color numbers 0-255
                    text.set_color((clade.color.red / 255, clade.color.green / 255, clade.color.blue / 255))

    def color_msa(seq):
        # Define nucleotide colors
        nucleotide_colors = {"A": "green", "C": "purple", "G": "orange", "T": "red", "-": "grey"}
        nucleotide_background_colors = {"A": "lightblue", "C": "lightblue", "G": "lightblue", "T": "lightblue", "-": "white"}
        color_vec = [nucleotide_colors[n.upper()] for n in seq]
        color_background_vec = [nucleotide_background_colors[n.upper()] for n in seq]

        return color_vec, color_background_vec

    if msa is not None:  # Display also msa sequences

        # Extract species names from the tree
        tree_species = {leaf.name for leaf in tree.get_terminals()}
        # Filter MSA blocks by tree species
#        num_tree_species = len(tree_species)
        # Get y-positions of leaves to align MSA correctly
        y_positions = get_tree_y_positions(tree)
        y_positions = {key.rstrip('0123456789'):y_positions[key] for key in y_positions}


        # Determine the longest sequence length across all records for consistent x-axis limits
        max_seq_length = max(len(record.seq) for block in msa for record in block)
        ax_msa.set_xlim(0, max_seq_length)
        ax_msa.set_ylim(-0.5, max(y_positions.values()) + 0.5)

        # Plot each MSA row
        found_ctr=0
        found_species=set()
        for block in msa:
            max_seq_length = max(len(record.seq) for record in block)  # make len specific for one block!!!
            msa_width = ax_msa.get_window_extent().width
            char_width = msa_width / max_seq_length
            display_msa = {species_name.rstrip('0123456789').replace("HL", ""): "-" * max_seq_length for species_name in tree_species}
            display_color = {species_name.rstrip('0123456789').replace("HL", ""): ['grey'] * max_seq_length for species_name in tree_species}
            display_background_color = {species_name.rstrip('0123456789').replace("HL", ""): ['white'] * max_seq_length for species_name in tree_species}


            for record in block:  # Fill matrix
                species_name = record.id.split('.')[0].rstrip('0123456789').replace("HL", "")
#                print("Looking for species=", species_name)
                if species_name in y_positions:
                    display_msa[species_name] = str(record.seq)
                    display_color[species_name], display_background_color[species_name] = color_msa(str(record.seq))
            break  # do only one block!

        # Now display:
        for genome_ver in tree_species:  # loop on species tree
            species_name = genome_ver.split('.')[0].rstrip('0123456789').replace("HL", "")
            y_pos = 1 - ((y_positions[species_name]+0.8) / (len(y_positions) + 0.8))  # inverts for top alignment

            # Calculate the total width of ax_msa
            # Plot each nucleotide as a rectangle with background color
            for j, nucleotide in enumerate(display_msa[species_name]):
                x_pos = j / max_seq_length #  * char_width   # Scale j to fit within plot width
                ax_msa.add_patch(patches.Rectangle((x_pos, y_pos - 0.5/(len(y_positions) + 0.8)), char_width, 1/(len(y_positions) + 0.8),
                                                   transform=ax_msa.transAxes, color=display_background_color[species_name][j]))  # , edgecolor='none'
                # Add nucleotide text
                ax_msa.text(x_pos + 0.5/max_seq_length, y_pos, nucleotide, ha='center', va='center', fontsize=10,
                    transform=ax_msa.transAxes, color=display_color[species_name][j])

#            tree_species_missing_in_msa = {t.rstrip('0123456789') for t in tree_species} - found_species
#            sequence = '-' * max_seq_length  # didn't find species in this region
#        print("Tree species missing in MSA for genomic region: ", tree_species_missing_in_msa)
#        print("Num. Species found in region=", found_ctr, " ; num in tree: ", len(tree_species))

        # Hide x-axis labels for MSA
        ax_msa.set_xticks([])
        ax_msa.set_yticks([])
        plt.subplots_adjust(wspace=0.05)  # Reduce spacing between subplots

    # Save the plot to a file
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

