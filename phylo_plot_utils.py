import matplotlib.pyplot as plt
from matplotlib import cm, colors
import numpy as np
from msa_utils import *
from matplotlib import gridspec
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec
import matplotlib.transforms as mtransforms
from Bio import Phylo, AlignIO

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
    y_leaf_positions = {leaf.name: idx for idx, leaf in enumerate(leaves)}
    inter = tree.get_nonterminals()
    y_internal_positions = {inter.name: idx+0.5 for idx, inter in enumerate(inter)}

    return y_leaf_positions, y_internal_positions



def color_tree(tree, subtree_name=None, values=None, msa=None, scores=None, good_positions=None, split_score='',
               cmap_name='rainbow', output_file="out_tree.png"):  # coolwarm
    """
    Generalized function to color a tree.
    - Either provide a subtree name to color that subtree.
    - Or provide a vector of numerical values to color branches with a heatmap.

    Parameters:
    - tree: A Bio.Phylo tree object.
    - output_file: Path to save the output image.
    - subtree_name: Name of the internal node whose subtree will be colored.
    - values: A dictionary or array of numerical values (one for each branch).
    - msa: A Bio.Align.MultipleSeqAlignment object containing the alignment (optional).
           If provided, it will be displayed with rows aligned to tree leaf nodes.
    - cmap_name: Name of the colormap to use for heatmap coloring.
    - output_file: Name of file to save the figure.
    """

    def draw_colored_rectangles(clade, x, y, color1, color2, width=0.02, height=0.02):
        """Draw two adjacent rectangles on a branch."""
        rect1 = patches.Rectangle((x, y), width / 2, height, linewidth=0, edgecolor=None, facecolor=color1)
        rect2 = patches.Rectangle((x + width / 2, y), width / 2, height, linewidth=0, edgecolor=None, facecolor=color2)
        ax_tree.add_patch(rect1)
        ax_tree.add_patch(rect2)

    def leaf_labels(clade):
        if clade.is_terminal():
            # Return the clade's name if terminal (leaf) and adjust color accordingly
            return clade.name
        return None

    # Prepare color map if using numerical values
    cmap = plt.get_cmap(cmap_name)
    # Get y-positions of leaves to align MSA correctly
    y_leaf_positions, y_inter_positions = get_tree_y_positions(tree)
    y_leaf_positions = {key.rstrip('0123456789'): y_leaf_positions[key] for key in y_leaf_positions}

    if msa is None:
        fig, ax_tree = plt.subplots(figsize=(10, 10))         # Draw the tree
    else:  # add msa axes
        if scores is None:
            # Set up the grid layout to place MSA on the right
            fig = plt.figure(figsize=(10, 8))
            gs = gridspec.GridSpec(1, 2, width_ratios=[4, 3])  # Adjust width ratios for desired spacing

            # Tree and colorbar axis on the left  ,  MSA plot axis on the right
            ax_tree, ax_msa = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])
        else:
            fig = plt.figure(figsize=(12, 10))  # Adjusted figure size
            gs = gridspec.GridSpec(2, 2, height_ratios=[9, 1], width_ratios=[4, 3])  # Adjusted GridSpec

            # Tree and colorbar axis on the left, MSA plot axis on the right
            ax_tree, ax_msa = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
            ax_tree.set_ylabel('')

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
        if isinstance(values, list):  # allow two values
 #           print("list values!")
            values_list = list(values[0].values())
#            print("vl:", values_list)
            values1 = values[0]
            values_keys = values[0].keys()
            values2 = values[1]
            values2_list = list(values[1].values())
 #           print("vl2:", values_list)
 #           print("all values:")
            all_values = values_list + values2_list
        else:
            values_list = list(values.values())
            values_keys = values.keys()
            all_values = values_list
            values1 = values


        print("all values again: ", all_values)
        # Normalize and get the colormap
        norm = colors.Normalize(vmin=np.min(all_values), vmax=np.max(all_values))
#        norm = colors.Normalize(vmin=np.min(values_list), vmax=np.max(values_list))
        cmap = plt.get_cmap(cmap_name)
        print("internal_positions: ", y_inter_positions)

        values_dict_simple = {k.name: values1[k] for k in values_keys}
        # Map each clade to its color based on the normalized values
        for clade in tree.get_nonterminals() + tree.get_terminals():
            if clade.name in values_dict_simple:
                color_value = norm(values_dict_simple[clade.name])
                branch_color = cmap(color_value)
                # Assign color to clade
                clade.color = colors.rgb2hex(branch_color)

#                if isinstance(values, list):  # allow two values, draw tree twice
#                    rate1 = values[0].get(clade, None)
#                    rate2 = values[1].get(clade, None)
#
#                    print("rates: ", rate1, rate2, " ; cn=", clade.name)
#                    if rate1 is not None and rate2 is not None:
#                        # Get colors for each rate value
#                        color1 = cmap(norm(rate1))
#                        color2 = cmap(norm(rate2))
#
#                        # Draw branch colored rectangles
#                        print("clade: ", clade, " ; clade branch-length:", clade.branch_length)
#                        print("clade y:", clade.y)
#                        x = clade.branch_length or 0
#                        if clade.name in y_inter_positions:
#                            y = y_inter_positions[clade.name] or 0  # Placeholder for y-position; adjust based on tree structure
#                        else:
#                            y = y_leaf_positions[clade.name.rstrip('0123456789')] or 0

                        # Plot the rectangles for each rate
#                        print("draw rectangle at ", x, y)
#                        draw_colored_rectangles(clade, x, y, color1, color2)
        if isinstance(values, list):  # allow two values, draw tree twice
            y_offset = 0.05
            offset_transform = mtransforms.Affine2D().translate(0, -y_offset) + ax_tree.transData

            # Draw the offset tree on the same axes with the applied transformation
#            Phylo.draw(offset_tree, do_show=False, axes=ax, label_func=None, branch_labels=None,
#                       branch_color_func=lambda c: c.color, transform=offset_transform)

            # Create a twin y-axis to overlay the second tree
            ax_tree_twin = ax_tree.twinx()
            min_y, max_y = ax_tree.get_ylim()
            ax_tree_twin.set_ylim(min_y - y_offset, max_y - y_offset)

            my_phylo_draw(tree, label_func=leaf_labels, values1=values1, values2=values2, y_offset=0.3, do_show=False, axes=ax_tree)

#            Phylo.draw(tree, label_func=leaf_labels, do_show=False, axes=ax_tree_twin)  # , transform=offset_transform)
            values_dict_simple = {k.name: values2[k] for k in values_keys}

            # Shift the y-axis limits to create space for the second tree

            offset_tree = copy.deepcopy(tree)
            for clade in offset_tree.get_nonterminals() + offset_tree.get_terminals():
                if clade.name in values_dict_simple:
                    color_value = norm(values_dict_simple[clade.name])
                    branch_color = cmap(color_value)
                    # Assign color to clade
                    clade.color = colors.rgb2hex(branch_color)
        #            clade.y += y_offset
            tree = copy.deepcopy(offset_tree)
            print("Draw second tree offset!! ")

    else:
        raise ValueError("Either subtree_name or values must be provided.")




    # Set the global title for the figure
    if not isinstance(split_score, str):
        split_score = round(split_score, 3)
    fig.suptitle("Best split score:"+str(split_score), fontsize=16)  # Show score

    if not isinstance(values, list):
        Phylo.draw(tree, label_func=leaf_labels, do_show=False, axes=ax_tree)
    ax_tree.set_ylabel('')

    # If coloring by numerical values, add a color bar
    if values is not None:
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])  # Set a dummy array for the colorbar
        if scores is None:
            cbar = plt.colorbar(sm, ax=ax_tree, location='left')
        else:
            ax_cbar = fig.add_subplot(gs[1, 0])  # Added subplot for horizontal colorbar
            cbar = plt.colorbar(sm, cax=ax_cbar, orientation='horizontal')  # Horizontal orientation
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
        if isinstance(msa, str):
            print("msa file: ", msa)
            msa = AlignIO.parse(msa, "maf")
        # Extract species names from the tree
        tree_species = {leaf.name for leaf in tree.get_terminals()}
        # Filter MSA blocks by tree species


        # Plot each MSA row
        for block in msa:
#            print("run block: ", block)
            max_seq_length = max(len(record.seq) for record in block)  # make len specific for one block!!!
            ax_msa.set_xlim(0, max_seq_length)
            ax_msa.set_ylim(-0.5, max(y_leaf_positions.values()) + 0.5)
            msa_width = ax_msa.get_window_extent().width
            char_width = msa_width / max_seq_length
            display_msa = {species_name.rstrip('0123456789').replace("HL", ""): "-" * max_seq_length for species_name in tree_species}
            display_color = {species_name.rstrip('0123456789').replace("HL", ""): ['grey'] * max_seq_length for species_name in tree_species}
            display_background_color = {species_name.rstrip('0123456789').replace("HL", ""): ['white'] * max_seq_length for species_name in tree_species}

            for record in block:  # Fill matrix
                species_name = record.id.split('.')[0].rstrip('0123456789').replace("HL", "")
#                print("Looking for species=", species_name)
                if species_name in y_leaf_positions:
                    display_msa[species_name] = str(record.seq)
                    display_color[species_name], display_background_color[species_name] = color_msa(str(record.seq))
            break  # do only one block!

        # Now display:
        for genome_ver in tree_species:  # loop on species tree
            species_name = genome_ver.split('.')[0].rstrip('0123456789').replace("HL", "")
            y_pos = 1 - ((y_leaf_positions[species_name]+0.8) / (len(y_leaf_positions) + 0.8))  # inverts for top alignment

            # Calculate the total width of ax_msa
            # Plot each nucleotide as a rectangle with background color
            for j, nucleotide in enumerate(display_msa[species_name]):
                x_pos = j / max_seq_length #  * char_width   # Scale j to fit within plot width
                ax_msa.add_patch(patches.Rectangle((x_pos, y_pos - 0.5/(len(y_leaf_positions) + 0.8)), char_width, 1/(len(y_leaf_positions) + 0.8),
                                                   transform=ax_msa.transAxes, color=display_background_color[species_name][j]))  # , edgecolor='none'
                # Add nucleotide text
                ax_msa.text(x_pos + 0.5/max_seq_length, y_pos, nucleotide, ha='center', va='center', fontsize=10,
                    transform=ax_msa.transAxes, color=display_color[species_name][j])

        # Hide x-axis labels for MSA
        ax_msa.set_xticks([])
        ax_msa.set_yticks([])
        plt.subplots_adjust(wspace=0.05)  # Reduce spacing between subplots

    if scores is not None:
        print("phylop_score_type=", type(scores), " ; len=", len(scores), "msa_len=", max_seq_length)
        ax_scores = fig.add_subplot(gs[1, 1])  # Added subplot for phyloP score plot
#        ax_scores.plot(scores, color='black')  # Plot with black line
        if good_positions is None:
            print("set range")
            good_positions = range(1, len(scores)+1)
        print("Good positions for plots: ", good_positions)
        ax_scores.bar(good_positions, scores, color='black')  # bar-plot
        ax_scores.set_xlim(0, len(scores) - 1)
        ax_scores.set_xlabel('Alignment Position')
        ax_scores.set_ylabel('Scores')
        ax_scores.set_yticks([])  # Simplified appearance by removing y-ticks


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


def my_phylo_draw_phydraw(
    tree,
    label_func=str,
    do_show=True,
    show_confidence=True,
    # For power users
    axes=None,
    branch_labels=None,
    label_colors=None,
    *args,
    **kwargs,
):
    """Plot the given tree using matplotlib (or pylab).
    The graphic is a rooted tree, drawn with roughly the same algorithm as
    draw_ascii.
    Additional keyword arguments passed into this function are used as pyplot
    options. The input format should be in the form of:
    pyplot_option_name=(tuple), pyplot_option_name=(tuple, dict), or
    pyplot_option_name=(dict).
    Example using the pyplot options 'axhspan' and 'axvline'::
        from Bio import Phylo, AlignIO
        from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
        constructor = DistanceTreeConstructor()
        aln = AlignIO.read(open('TreeConstruction/msa.phy'), 'phylip')
        calculator = DistanceCalculator('identity')
        dm = calculator.get_distance(aln)
        tree = constructor.upgma(dm)
        Phylo.draw(tree, axhspan=((0.25, 7.75), {'facecolor':'0.5'}),
        ... axvline={'x':0, 'ymin':0, 'ymax':1})
    Visual aspects of the plot can also be modified using pyplot's own functions
    and objects (via pylab or matplotlib). In particular, the pyplot.rcParams
    object can be used to scale the font size (rcParams["font.size"]) and line
    width (rcParams["lines.linewidth"]).
    :Parameters:
        label_func : callable
            A function to extract a label from a node. By default this is str(),
            but you can use a different function to select another string
            associated with each node. If this function returns None for a node,
            no label will be shown for that node.
        do_show : bool
            Whether to show() the plot automatically.
        show_confidence : bool
            Whether to display confidence values, if present on the tree.
        axes : matplotlib/pylab axes
            If a valid matplotlib.axes.Axes instance, the phylogram is plotted
            in that Axes. By default (None), a new figure is created.
        branch_labels : dict or callable
            A mapping of each clade to the label that will be shown along the
            branch leading to it. By default this is the confidence value(s) of
            the clade, taken from the ``confidence`` attribute, and can be
            easily toggled off with this function's ``show_confidence`` option.
            But if you would like to alter the formatting of confidence values,
            or label the branches with something other than confidence, then use
            this option.
        label_colors : dict or callable
            A function or a dictionary specifying the color of the tip label.
            If the tip label can't be found in the dict or label_colors is
            None, the label will be shown in black.
    """
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        try:
            import pylab as plt
        except ImportError:
            raise ValueError(
                "Install matplotlib or pylab if you want to use draw."
            ) from None
    import matplotlib.collections as mpcollections
    # Arrays that store lines for the plot of clades
    horizontal_linecollections = []
    vertical_linecollections = []
    # Options for displaying branch labels / confidence
    def conf2str(conf):
        if int(conf) == conf:
            return str(int(conf))
        return str(conf)
    if not branch_labels:
        if show_confidence:
            def format_branch_label(clade):
                try:
                    confidences = clade.confidences
                    # phyloXML supports multiple confidences
                except AttributeError:
                    pass
                else:
                    return "/".join(conf2str(cnf.value) for cnf in confidences)
                if clade.confidence is not None:
                    return conf2str(clade.confidence)
                return None
        else:
            def format_branch_label(clade):
                return None
    elif isinstance(branch_labels, dict):
        def format_branch_label(clade):
            return branch_labels.get(clade)
    else:
        if not callable(branch_labels):
            raise TypeError(
                "branch_labels must be either a dict or a callable (function)"
            )
        format_branch_label = branch_labels
    # options for displaying label colors.
    if label_colors:
        if callable(label_colors):
            def get_label_color(label):
                return label_colors(label)
        else:
            # label_colors is presumed to be a dict
            def get_label_color(label):
                return label_colors.get(label, "black")
    else:
        def get_label_color(label):
            # if label_colors is not specified, use black
            return "black"
    # Layout
    def get_x_positions(tree):
        """Create a mapping of each clade to its horizontal position.
        Dict of {clade: x-coord}
        """
        depths = tree.depths()
        # If there are no branch lengths, assume unit branch lengths
        if not max(depths.values()):
            depths = tree.depths(unit_branch_lengths=True)
        return depths
    def get_y_positions(tree):
        """Create a mapping of each clade to its vertical position.
        Dict of {clade: y-coord}.
        Coordinates are negative, and integers for tips.
        """
        maxheight = tree.count_terminals()
        # Rows are defined by the tips
        heights = {
            tip: maxheight - i for i, tip in enumerate(reversed(tree.get_terminals()))
        }
        # Internal nodes: place at midpoint of children
        def calc_row(clade):
            for subclade in clade:
                if subclade not in heights:
                    calc_row(subclade)
            # Closure over heights
            heights[clade] = (
                heights[clade.clades[0]] + heights[clade.clades[-1]]
            ) / 2.0
        if tree.root.clades:
            calc_row(tree.root)
        return heights
    x_posns = get_x_positions(tree)
    y_posns = get_y_positions(tree)
    # The function draw_clade closes over the axes object
    if axes is None:
        fig = plt.figure()
        axes = fig.add_subplot(1, 1, 1)
    elif not isinstance(axes, plt.matplotlib.axes.Axes):
        raise ValueError(f"Invalid argument for axes: {axes}")
    def draw_clade_lines(
        use_linecollection=False,
        orientation="horizontal",
        y_here=0,
        x_start=0,
        x_here=0,
        y_bot=0,
        y_top=0,
        color="black",
        lw=".1",
    ):
        """Create a line with or without a line collection object.
        Graphical formatting of the lines representing clades in the plot can be
        customized by altering this function.
        """
        if not use_linecollection and orientation == "horizontal":
            axes.hlines(y_here, x_start, x_here, color=color, lw=lw)
        elif use_linecollection and orientation == "horizontal":
            horizontal_linecollections.append(
                mpcollections.LineCollection(
                    [[(x_start, y_here), (x_here, y_here)]], color=color, lw=lw
                )
            )
        elif not use_linecollection and orientation == "vertical":
            axes.vlines(x_here, y_bot, y_top, color=color)
        elif use_linecollection and orientation == "vertical":
            vertical_linecollections.append(
                mpcollections.LineCollection(
                    [[(x_here, y_bot), (x_here, y_top)]], color=color, lw=lw
                )
            )
    def draw_clade(clade, x_start, color, lw):
        """Recursively draw a tree, down from the given clade."""
        x_here = x_posns[clade]
        y_here = y_posns[clade]
        # phyloXML-only graphics annotations
        if hasattr(clade, "color") and clade.color is not None:
            color = clade.color.to_hex()
        if hasattr(clade, "width") and clade.width is not None:
            lw = clade.width * plt.rcParams["lines.linewidth"]
        # Draw a horizontal line from start to here
        draw_clade_lines(
            use_linecollection=True,
            orientation="horizontal",
            y_here=y_here,
            x_start=x_start,
            x_here=x_here,
            color=color,
            lw=lw,
        )
        # Add node/taxon labels
        label = label_func(clade)
        if label not in (None, clade.__class__.__name__):
            axes.text(
                x_here,
                y_here,
                f" {label}",
                verticalalignment="center",
                color=get_label_color(label),
            )
        # Add label above the branch (optional)
        conf_label = format_branch_label(clade)
        if conf_label:
            axes.text(
                0.5 * (x_start + x_here),
                y_here,
                conf_label,
                fontsize="small",
                horizontalalignment="center",
            )
        if clade.clades:
            # Draw a vertical line connecting all children
            y_top = y_posns[clade.clades[0]]
            y_bot = y_posns[clade.clades[-1]]
            # Only apply widths to horizontal lines, like Archaeopteryx
            draw_clade_lines(
                use_linecollection=True,
                orientation="vertical",
                x_here=x_here,
                y_bot=y_bot,
                y_top=y_top,
                color=color,
                lw=lw,
            )
            # Draw descendents
            for child in clade:
                draw_clade(child, x_here, color, lw)
    draw_clade(tree.root, 0, "k", plt.rcParams["lines.linewidth"])
    # If line collections were used to create clade lines, here they are added
    # to the pyplot plot.
    for i in horizontal_linecollections:
        axes.add_collection(i)
    for i in vertical_linecollections:
        axes.add_collection(i)
    # Aesthetics
    try:
        name = tree.name
    except AttributeError:
        pass
    else:
        if name:
            axes.set_title(name)
    axes.set_xlabel("branch length")
    axes.set_ylabel("taxa")
    # Add margins around the tree to prevent overlapping the axes
    xmax = max(x_posns.values())
    axes.set_xlim(-0.05 * xmax, 1.25 * xmax)
    # Also invert the y-axis (origin at the top)
    # Add a small vertical margin, but avoid including 0 and N+1 on the y axis
    axes.set_ylim(max(y_posns.values()) + 0.8, 0.2)
    # Parse and process key word arguments as pyplot options
    for key, value in kwargs.items():
        try:
            # Check that the pyplot option input is iterable, as required
            list(value)
        except TypeError:
            raise ValueError(
                'Keyword argument "%s=%s" is not in the format '
                "pyplot_option_name=(tuple), pyplot_option_name=(tuple, dict),"
                " or pyplot_option_name=(dict) " % (key, value)
            ) from None
        if isinstance(value, dict):
            getattr(plt, str(key))(**dict(value))
        elif not (isinstance(value[0], tuple)):
            getattr(plt, str(key))(*value)
        elif isinstance(value[0], tuple):
            getattr(plt, str(key))(*value[0], **dict(value[1]))
    if do_show:
        plt.show()


def my_phylo_draw(
        tree,
        values1,
        values2,
        label_func=str,
        do_show=True,
        show_confidence=True,
        # For power users
        axes=None,
        branch_labels=None,
        label_colors=None,
        y_offset=0.5,  # Vertical offset for the second tree
        cmap_name='rainbow',  # Colormap for coloring branches
        *args,
        **kwargs,
):
    """Plot the given tree using matplotlib with two color schemes for branch rates.

    Parameters:
        values1 : dict
            A dictionary mapping clade names to rate values for the first color scheme.
        values2 : dict
            A dictionary mapping clade names to rate values for the second color scheme.
        y_offset : float
            Vertical offset for the second tree to separate rate vectors visually.
    """
    import matplotlib.pyplot as plt
    import matplotlib.collections as mpcollections
    from matplotlib import colors

    # Get colormap and normalization based on combined rate vectors
    cmap = plt.get_cmap(cmap_name)
    all_values = list(values1.values()) + list(values2.values())
    norm = colors.Normalize(vmin=min(all_values), vmax=max(all_values))

    # Arrays that store lines for the plot of clades
    horizontal_linecollections = []
    vertical_linecollections = []

    # Custom color and label formatting functions
    def format_branch_label(clade):
        return str(clade.confidence) if show_confidence and clade.confidence else None

    def get_color_for_clade(clade, values):
        """Get color for a clade based on its rate value."""
        if clade.name in values:
            rate = values[clade.name]
            return cmap(norm(rate))
        return "black"  # Default color if clade not in values

    # Layout
    def get_x_positions(tree):
        depths = tree.depths(unit_branch_lengths=not max(tree.depths().values()))
        return depths

    def get_y_positions(tree):
        maxheight = tree.count_terminals()
        heights = {tip: maxheight - i for i, tip in enumerate(reversed(tree.get_terminals()))}

        def calc_row(clade):
            for subclade in clade:
                if subclade not in heights:
                    calc_row(subclade)
            heights[clade] = (heights[clade.clades[0]] + heights[clade.clades[-1]]) / 2.0

        if tree.root.clades:
            calc_row(tree.root)
        return heights

    x_posns = get_x_positions(tree)
    y_posns = get_y_positions(tree)

    # Setup Matplotlib axes
    if axes is None:
        fig, axes = plt.subplots(figsize=(10, 8))

    # Function to draw a tree
    def draw_tree(clade, x_start, values, y_adjust=0):
        """Recursively draw branches for a tree with given values and y-offset."""
        x_here = x_posns[clade]
        y_here = y_posns[clade] + y_adjust
        color = get_color_for_clade(clade, values)

        # Draw horizontal branch line
        axes.plot([x_start, x_here], [y_here, y_here], color=color, lw=2)

        # Add label at the branch
        label = label_func(clade)
        if label:
            axes.text(x_here, y_here, f" {label}", verticalalignment="center", color=color)

        # Draw vertical line connecting children
        if clade.clades:
            y_top = y_posns[clade.clades[0]] + y_adjust
            y_bot = y_posns[clade.clades[-1]] + y_adjust
            axes.plot([x_here, x_here], [y_bot, y_top], color=color, lw=2)

            # Recursively draw each child
            for child in clade:
                draw_tree(child, x_here, values, y_adjust)

    # Draw the first tree using `values1`
    draw_tree(tree.root, 0, values1, y_adjust=0)

    # Draw the second tree using `values2` with a slight vertical offset
    draw_tree(tree.root, 0, values2, y_adjust=-y_offset)

    # Color bar to represent the values scale
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=axes, orientation="horizontal", pad=0.05)
    cbar.set_label("Rate")

    # Aesthetics
    axes.set_xlabel("Branch Length")
    axes.set_ylabel("Taxa")
    axes.set_ylim(max(y_posns.values()) + 1, min(y_posns.values()) - y_offset - 1)
    xmax = max(x_posns.values())
    axes.set_xlim(-0.05 * xmax, 1.05 * xmax)

    if do_show:
        plt.show()