import matplotlib.pyplot as plt
from matplotlib import cm, colors
from Bio import Phylo
import numpy as np


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


def color_tree(tree, subtree_name=None, values=None, msa=None, cmap_name='viridis', output_file="out_tree.png"):
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
        # Normalize the values to [0,1] for the color map
        norm = colors.Normalize(vmin=np.min(values), vmax=np.max(values))
        branch_colors = cmap(norm(values))

        # Assign colors to each branch based on the values
        for clade, branch_color in zip(tree.get_nonterminals() + tree.get_terminals(), branch_colors):
            # Convert the RGB tuple (from colormap) to a hex string for matplotlib
            clade.color = colors.rgb2hex(branch_color)

    else:
        raise ValueError("Either subtree_name or values must be provided.")

        # Custom label function to return the name for terminal nodes, and set leaf color

    def leaf_labels(clade):
        if clade.is_terminal():
            # Return the clade's name if terminal (leaf) and adjust color accordingly
            return clade.name
        return None

    # Draw the tree
    fig = plt.figure(figsize=(10, 10))
    axes = fig.add_subplot(1, 1, 1)

    # Draw the tree with colored branches
    Phylo.draw(tree, axes=axes, label_func=leaf_labels,  do_show=False)

    # If coloring by numerical values, add a color bar
    if values is not None:
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])  # Set a dummy array for the colorbar
        cbar = plt.colorbar(sm, ax=axes)
        cbar.set_label('Branch Values')
#    else:
    for text in axes.texts:  # axes.texts contains all the text objects
        leaf_name = text.get_text()
        for clade in tree.get_terminals():
            if clade.name == leaf_name.strip():
#                print("IF: Setting color: ", leaf_name , ", " , clade.color, type(clade.color))
                # Set the color of the text based on the clade's assigned color
                if isinstance(clade.color, str):
                    text.set_color(clade.color)
                else:
                    text.set_color((clade.color.red / 255, clade.color.green / 255, clade.color.blue / 255))

    # Save the plot to a file
  #  print("SAVING!!!")
    plt.savefig(output_file)
    plt.close()


def color_tree_with_msa(tree, values_vector, msa=None):
    """
    Draws a phylogenetic tree with colored branches according to `values_vector`
    and displays a multiple sequence alignment (MSA) next to the tree if provided.

    Parameters:
    - tree: A Bio.Phylo tree object to be drawn.
    - values_vector: A list of color values corresponding to each branch in the tree.
    - msa: A Bio.Align.MultipleSeqAlignment object containing the alignment (optional).
           If provided, it will be displayed with rows aligned to tree leaf nodes.
    """
    # Create a color mapping for the branches based on `values_vector`
    unique_values = sorted(set(values_vector))
    colors = plt.cm.tab20(np.linspace(0, 1, len(unique_values)))
    color_map = {value: colors[i] for i, value in enumerate(unique_values)}

    # Draw the phylogenetic tree
    fig, ax_tree = plt.subplots(figsize=(10, 6))
    ax_tree.set_title("Phylogenetic Tree with Colored Branches")

    # Assign colors to branches
    for clade, color_value in zip(tree.find_clades(), values_vector):
        clade.color = color_map[color_value]

    # Draw the tree
    Phylo.draw(tree, axes=ax_tree, branch_labels=None, do_show=False)

    # If an MSA is provided, draw it alongside the tree
    if msa:
        # Collect leaf names in the order they appear in the tree drawing
        leaves = tree.get_terminals()
        leaf_names = [leaf.name for leaf in leaves]

        # Create a dictionary mapping leaf names to their MSA sequences
        msa_dict = {record.id: str(record.seq) for record in msa}

        # Create an alignment matrix to display the sequences alongside the tree
        msa_matrix = [list(msa_dict[leaf_name]) if leaf_name in msa_dict else [''] * msa.get_alignment_length() for
                      leaf_name in leaf_names]

        # Convert MSA matrix to numpy for easier display
        msa_array = np.array(msa_matrix)

        # Create a new axes for the MSA display
        ax_msa = fig.add_axes([0.8, 0.1, 0.15, 0.8])  # Adjust position and size as needed
        ax_msa.set_title("MSA")

        # Display the MSA matrix as an image
        # Color mapping for nucleotides or amino acids
        msa_cmap = plt.cm.Pastel1
        ax_msa.imshow(msa_array == 'A', cmap=msa_cmap, aspect='auto',
                      interpolation='none')  # Example with 'A'; update to generalize

        # Adjust the axis
        ax_msa.set_xticks(range(msa.get_alignment_length()))
        ax_msa.set_yticks(range(len(leaf_names)))
        ax_msa.set_yticklabels(leaf_names, fontsize=8)
        ax_msa.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)

    plt.show()


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

        for clade in tree.find_clades():
            # Assign the color only if the clade is a terminal node and in the current subtree
            if clade.is_terminal() and clade.name in subtree_species:
                values_dict[clade] = color_index

    # Assign color 0 to all branches not part of any subtree
    for clade in tree.find_clades():
        if clade not in values_dict:
            values_dict[clade] = 0

    # Convert values_dict to a list vector for the color_tree function
    values_vector = [values_dict[clade] for clade in tree.find_clades()]

    return values_vector

