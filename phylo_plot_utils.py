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


def color_tree(tree, subtree_name=None, values=None, cmap_name='viridis', output_file="out_tree.png"):
    """
    Generalized function to color a tree.
    - Either provide a subtree name to color that subtree.
    - Or provide a vector of numerical values to color branches with a heatmap.

    Parameters:
    - tree: A Bio.Phylo tree object.
    - output_file: Path to save the output image.
    - subtree_name: Name of the internal node whose subtree will be colored.
    - values: A list or array of numerical values (one for each branch).
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


# Color a tree with subtree in different color
def color_subtree(tree, internal_node_name, output_file):
    # Function to color branches of a clade and its subclades
    def color_clades(clade, color):
        for subclade in clade.clades:
            subclade.color = color  # Add color to each clade
            color_clades(subclade, color)  # Recursive coloring

    # Check if internal node exists in tree
    found_clade = None
    for clade in tree.find_clades():
        if clade.name == internal_node_name:
            found_clade = clade
            break

    if found_clade is None:
        raise ValueError(f"Internal node {internal_node_name} not found in the tree.")

    # Color the entire tree (default color)
    for clade in tree.get_terminals() + tree.get_nonterminals():
        clade.color = 'black'  # Set the default color to black for all nodes

    # Color the subtree of the found internal node
    color_clades(found_clade, 'red')  # Color the subtree in red

    # Custom label function to return the name for terminal nodes, and set leaf color
    def leaf_labels(clade):
        if clade.is_terminal():
            # Return the clade's name if terminal (leaf) and adjust color accordingly
            return clade.name
        return None

    # Drawing the tree with label_func and color adjustments
    fig = plt.figure(figsize=(10, 10))
    axes = fig.add_subplot(1, 1, 1)

    # Draw the tree
    Phylo.draw(tree, axes=axes, label_func=leaf_labels, do_show=False)

    # Now adjust the text colors of the leaves after drawing
    for text in axes.texts:  # axes.texts contains all the text objects
        leaf_name = text.get_text()
     #   print("text: ", text, " leaf: ", leaf_name)
     #   print("Terminals: ")
     #   print(tree.get_terminals())
        for clade in tree.get_terminals():
     #       print("clade name: ", clade, " ; ", clade.name, " ; leaf_name: ", leaf_name, "..", leaf_name[1:], " equal ? ", \
     #             clade.name == leaf_name[1:], " lens: ", len(str(clade)), len(clade.name), " ; ", len(leaf_name))

            if clade.name == leaf_name.strip():
      #          print("IF: Setting color: ", leaf_name , ", " , clade.color, type(clade.color))
                # Set the color of the text based on the clade's assigned color
                if isinstance(clade.color, str):
                    text.set_color(clade.color)
                else:
    #                print("Tuple: ", tuple(clade.color))
    #                norm_color = [x/256 for x in tuple(clade.color)]
                    text.set_color((clade.color.red / 255, clade.color.green / 255, clade.color.blue / 255))



    # Save the plot to a file
    plt.savefig(output_file)
    plt.close()
