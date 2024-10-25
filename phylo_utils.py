import subprocess
from Bio import Phylo
# from phasty import *
from phylo_plot_utils import *
import numpy as np
import os
import copy
import importlib


def split_tree_at_root(tree):
    # Get the children of the root
    root = tree.root
    children = root.clades  # Get the two child clades (subtrees)

    if len(children) != 2:
        raise ValueError("Tree is not bifurcating at the root, cannot split into two subtrees.")

    # Create two subtrees based on the root's children
    subtree_1, subtree_2 = children

    # Rename the internal nodes of the new subtrees
    subtree_1.name, subtree_2.name = "Subtree_1", "Subtree_2"

    # Create two separate trees with each subtree
    tree_1 = Phylo.BaseTree.Tree(root=subtree_1)
    tree_2 = Phylo.BaseTree.Tree(root=subtree_2)

    return tree_1, tree_2


def extract_and_prune_model(input_file, species_to_remove, output_file):
    """
    Extract the Newick tree from the model, prune the specified species,
    and replace the tree back in the model.

    Parameters:
    - input_file: Path to the input model file containing substitution model and tree.
    - species_to_remove: A list of species names to remove from the tree.
    - output_file: Path to save the output model with the pruned tree.
    """
    # Step 1: Read the model and extract the tree
    with open(input_file, "r") as infile:
        lines = infile.readlines()

    # Find the line that starts with "TREE:" and extract the Newick tree
    tree_line_index = None
    for idx, line in enumerate(lines):
        if line.startswith("TREE:"):
            tree_line_index = idx
            break

    if tree_line_index is None:
        raise ValueError("No TREE line found in the model file.")

    # Extract the Newick tree (after "TREE: ")
    newick_tree = lines[tree_line_index].strip().split("TREE:")[1].strip()
    # Step 2: Load the Newick tree into Biopython and prune the species

    tree = Phylo.read(input_file, "newick")
#    tree = Phylo.read(newick_tree, "newick")

    for species in species_to_remove:
        clade_to_remove = tree.find_any(name=species)
        if clade_to_remove:
            tree.prune(clade_to_remove)
    # Step 3: Write the pruned tree back into the original model
    # Convert the pruned tree back to a Newick string
    pruned_tree_newick = tree.format("newick").strip()

    # Replace the old tree line with the new pruned tree
    lines[tree_line_index] = f"TREE: {pruned_tree_newick}\n"

    # Step 4: Write the updated model back to the output file
    with open(output_file, "w") as outfile:
        outfile.writelines(lines)


def prune_tree(tree, species_to_remove):
    """
    Prune the tree by removing the species listed in species_to_remove.

    Parameters:
    - tree: A Bio.Phylo tree object.
    - species_to_remove: A list of species names to remove from the tree.

    Returns:
    - A pruned tree object with the specified species removed.
    """
    for species in species_to_remove:
        # Find and prune each species from the tree
        clade_to_remove = tree.find_any(name=species)
        if clade_to_remove:
            tree.prune(clade_to_remove)

    return tree


def get_induced_subtree(tree, root_name):
    """
    Returns the subtree of `tree` with `root_name` as the root.

    Parameters:
    - tree: A Bio.Phylo tree object containing all species.
    - root_name: The name of the node to be used as the root of the induced subtree.

    Returns:
    - A new Bio.Phylo tree object representing the induced subtree with `root_name` as the root.
    """
    # Clone the original tree to avoid modifying it
    subtree = copy.deepcopy(tree)

    # Find the specified node in the tree
    root_clade = subtree.find_any(name=root_name)
    if root_clade is None:
        raise ValueError(f"No clade with the name '{root_name}' found in the tree.")

    # Trim the tree down to the specified root node and its descendants
    subtree.root = root_clade
    subtree.rooted = True

    # Detach all nodes above the specified root to make it the new root
    subtree.collapse_all(lambda clade: clade is not root_clade)

    return subtree


def get_complementary_tree(tree, sub_tree):
    """
    Returns the complementary tree of `tree` by removing all species in `sub_tree`.

    Parameters:
    - tree: A Bio.Phylo tree object containing all species.
    - sub_tree: A Bio.Phylo tree object containing a subset of species in `tree`.

    Returns:
    - A new tree object representing the complementary tree.
    """
    if isinstance(sub_tree, str):
        sub_tree = get_induced_subtree(tree, sub_tree)
    print("sub_tree is: ", sub_tree)
    print("terminals are: ", sub_tree.get_terminals())
    # Get the names of all species (terminal nodes) in sub_tree
    sub_tree_species = {leaf.name for leaf in sub_tree.get_terminals()}

    # Clone the original tree to avoid modifying it directly
    complementary_tree = copy.deepcopy(tree)

    # Prune each species in sub_tree from complementary_tree
    for species in sub_tree_species:
        clade_to_remove = complementary_tree.find_any(name=species)
        if clade_to_remove:
            complementary_tree.prune(clade_to_remove)

    return complementary_tree


def read_phylop_output_fixed_step(output_file, plot_flag=False):
    """
    Reads PhyloP output in fixedStep format into a list of dictionaries.

    Parameters:
    - output_file: Path to the PhyloP output file.

    Returns:
    - A list of dictionaries, each containing the chromosome, position, and test statistics.
    """
    phylop_scores = []
    chrom, start, position = [None]*3
    step = 1

    with open(output_file, 'r') as f:
        for line in f:
            # Skip comment lines or headers
            if line.startswith("#") or not line.strip():
                continue

            # Parse the fixedStep line to get chrom and start position
            if line.startswith("fixedStep"):
                # Example: fixedStep chrom=chr22_GL383583v2_alt start=1 step=1
                parts = line.strip().split()
                chrom = parts[1].split('=')[1]  # Extract chromosome
                start = int(parts[2].split('=')[1])  # Extract start position
                step = int(parts[3].split('=')[1])  # Extract step size
                position = start  # Initialize position
                continue

            # Process the data lines
            # Example data line format: 0.00000 -0.43885 -0.00000 0.00000 1.00000
            values = line.strip().split()

            # We expect five columns: scale, deriv, subderiv, teststat, pval
            if len(values) == 5:
                scale, deriv, subderiv, teststat, pval = map(float, values)

                # Append the current record to the list
                phylop_scores.append({
                    "chromosome": chrom,
                    "position": position,
                    "scale": scale,
                    "deriv": deriv,
                    "subderiv": subderiv,
                    "teststat": teststat,
                    "pval": pval
                })

                # Increment the position by step
                position += step

    if plot_flag:
        plot_phylop_scores(phylop_scores, output_file.replace(".out", ".png"))

    return phylop_scores

# Read output file into a list of  dictionaries
def read_phylop_output(output_file, plot_flag = False):
    """
    Reads PhyloP base-by-base output into a list of dictionaries.

    Parameters:
    - output_file: Path to the PhyloP output file.

    Returns:
    - A list of dictionaries, each containing the chromosome, position, and score.
    """
    phylop_scores = []

    with open(output_file, 'r') as f:
        for line in f:
            # Skip comment lines or headers
            if line.startswith("#") or not line.strip():
                continue
#
#            # Split the line into components (assumes format: chr, pos, score)
            print("line is: ", line)
            cur_line = line.strip().split()
            print(cur_line)
            chrom, position, score = line.strip().split()
#
#            # Convert position and score to appropriate types
            position = int(position)
            score = float(score)
#
#            # Store each result as a dictionary
            phylop_scores.append({
                "chromosome": chrom,
                "position": position,
                "score": score
            })

    if plot_flag:
        plot_phylop_scores(phylop_scores, output_file.replace(".out", ".png"))

    return phylop_scores


# Run phylop using the command line
def run_phylop_linux(tree_file, msa_file, output_file, sub_tree="", method="SCORE", mode="CONACC", positions_file=None, read_output=False):
    """
    Run PhyloP from PHAST installed on WSL (Linux subsystem) using Python's subprocess module.

    Parameters:
    - tree_file: Path to the tree file in Newick format (path within WSL).
    - msa_file: Path to the multiple sequence alignment file (path within WSL).
    - output_file: Path to save the PhyloP output (path within WSL).
    - mode: PhyloP mode ('CONACC', 'CON', or 'ACC'). Default is 'CONACC'.
    """

    tree_file = "'" + tree_file + "'"
    msa_file = "'" + msa_file + "'"
    output_file = "'" + output_file + "'"
    if positions_file is None:
        pos_str = ""
    else:
        positions_file = "'" + positions_file + "'"
        pos_str = f"--only-sites {positions_file} "
    command = f"/mnt/c/Code/Github/phast/bin/phyloP --method {method} --subtree {sub_tree} --mode {mode} --base-by-base " + pos_str + \
              f" {tree_file}  {msa_file} > {output_file}"

    print("Run phylop command: ")
    print(command)
    # Call the WSL command from Windows
    subprocess.run(command, shell=True)  # now wsl !!!
#    with open(output_file, "w") as outfile:
#        subprocess.run(["wsl", command], stdout=outfile, shell=True)

    if read_output:
        print("Read phylop output: ", output_file, output_file.replace("'", ""))
        return read_phylop_output_fixed_step(output_file.replace("'", ""))  # read with windows format
    else:
        return None


# New: split tree recursively and output subtree with different rates !
def fit_tree_rates(tree_file, msa_file, output_file, fdr_alpha = 0.1,  method="SCORE", mode="CONACC", plot_tree=False):
    tree = Phylo.read(tree_file, "newick")
    root = tree.root
    children = root.clades  # Get the two child clades (subtrees)
#    phylop_scores = run_phylop_linux(tree_file, msa_file, output_file, sub_tree=children[0], method="SCORE", mode="CONACC", read_output=True)

    # Decide if significant or not
    output_dir = os.path.dirname(output_file)
    best_subtree, best_score = find_best_rate_split_in_tree(tree_file, msa_file, output_dir, "greedy")
    print("Best Subtree: ", best_subtree)
    print(type(best_subtree))
    print("Best score: ", best_score)
    comp_best_tree = get_complementary_tree(tree, best_subtree) # Fit rates for the two trees
    if plot_tree:
        print("Plotting! Saving in: ", output_file.replace(".out", ".png"))
        subtree_rates = subtrees_to_color_vector(tree, [get_induced_subtree(tree, best_subtree)])  # Compute rates and color
        color_tree(tree, values=subtree_rates, cmap_name='viridis', output_file=output_file.replace(".out", ".png"))

    return best_subtree, best_score


# Split the tree such that the contrast in rate between the two halves is maximal.
# Default: Use greedy approach with binary search. Alternative is to enumerate all possible splits
def find_best_rate_split_in_tree(tree_file, msa_file, output_dir, method="binary"):
    tree = Phylo.read(tree_file,  "newick")
    print("Start best rate split, output-dir=", output_dir)
    if method == "greedy":  # enumerate over all possible splits

        best_subtree = None
        best_score = float('-inf')

        # Enumerate all internal nodes (subtrees) AND leaf nodes (species)
        for clade in tree.find_clades():
            # Get the name of the clade (species or internal node)
            subtree_name = clade.name or f"clade_{clade.confidence}"  # Use name or confidence if unnamed

            # Skip if the subtree has no name (rare but possible)
            if not subtree_name:
                continue

            # Create a temporary output file for PhyloP results
            print("output dir: ", output_dir)
            output_file = os.path.join(output_dir, f"phylop_{subtree_name}.out")

            print("Run phylop outputfile=", output_file)
            # Run PhyloP for the current subtree (species or subtree)
            phylop_scores = run_phylop_linux(tree_file, msa_file, output_file, \
                                             sub_tree=subtree_name, mode="CONACC", method="SCORE", read_output=True)
            print("Finished phylop outputfile=", output_file)
            # Process the PhyloP output and calculate the score difference
            cur_score = sum([ x['teststat'] for x in phylop_scores])

            # Track the subtree/species with the maximum score
            if cur_score > best_score:
                best_score = cur_score
                best_subtree = subtree_name
                print("Found better!!!")
                print(best_subtree, best_score)

        return best_subtree, best_score


