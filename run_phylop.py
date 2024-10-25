# Run phylop routines
from Bio import Phylo
import matplotlib.pyplot as plt
import os

import numpy as np
from matplotlib import cm, colors
from phylo_utils import *
from sys import platform

# Operations to do:
plot_trees = False
run_phylop = False
run_phylop_recursive = True


# Figure out operating system
if platform == "linux" or platform == "linux2":
    run_os = "linux"     # linux
    g_dir = "/mnt/g/My Drive/"
elif platform == "darwin":
    run_os = "darwin"
    # OS X
elif platform == "win32" or platform == "win64":
    run_os = "windows" # Windows...
    g_dir = "G://My Drive/"



# Set main data directory:
data_dir = g_dir + "Data/hg38_470Vert"
tree_file = data_dir + "/phylop30/hg38.phyloP30way_named.mod"
msa_file = data_dir + "/multiz30/chr22_GL383583v2_alt.maf"
output_dir = data_dir + "/output"

# Read the tree from a Newick file
tree = Phylo.read(tree_file, "newick")
# Prune tree and save
missing_species_in_msa = ['ponAbe2', 'papAnu3', 'nasLar1', 'rhiRox1', 'rhiBie1', 'tarSyr2', 'micMur3', 'proCoq1', 'eulMac1', 'eulFla1', 'mm10', 'dasNov3']
# pruned_tree = prune_tree(tree, missing_species_in_msa)
pruned_tree_file = tree_file[:-4] + "_pruned.mod"
# Phylo.write(pruned_tree, pruned_tree_file, "newick")
print("Prune tree from mod file:")
extract_and_prune_model(tree_file, missing_species_in_msa, pruned_tree_file)
print("READ Prune tree from mod file:")
pruned_tree = Phylo.read(pruned_tree_file, "newick")


if plot_trees:
    # Print the entire tree structure
    Phylo.draw_ascii(tree)

    # Access internal nodes and leaves
    for clade in tree.find_clades():
        if clade.is_terminal():
            print(f"Leaf: {clade.name}")
        else:
            print(f"Internal Node: {clade}")
    for clade in tree.find_clades():
        print(f"Node: {clade}, Branch length: {clade.branch_length}")

    # Example usage:
    # tree = Phylo.read("your_tree_file.newick", "newick")

    print("Color trees:")
    color_tree(tree, subtree_name="hg38-calJac3", output_file=data_dir + "/output/trees/output_tree_hg38-calJac3_subtree.png")
    color_tree(tree, subtree_name="calJac3-saiBol1", output_file=data_dir + "/output/trees/output_tree_calJac3-saiBol1_subtree.png")
    color_tree(pruned_tree, subtree_name="hg38-calJac3", output_file=data_dir + "/output/trees/output_pruned_tree_hg38-calJac3_subtree.png")
    color_tree(pruned_tree, subtree_name="hg38-rheMac8", output_file=data_dir + "/output/trees/output_pruned_tree_hg38-rheMac8_subtree.png")

    values = np.random.random(len(tree.get_terminals()) + len(tree.get_nonterminals()))
    color_tree(tree, values=values, output_file=data_dir + "/output/trees/output_tree_hg38_randcolors.png")


# Run phylop:
# Example usage:
# Paths must be within the WSL environment (e.g., /home/username/data/)
# tree_file = "/mnt/c/path/to/tree_file.newick"
# msa_file = "/mnt/c/path/to/msa_file.maf"
sub_tree = "calJac3-saiBol1"
output_phylop_file = data_dir + "/output/phylop/phylop_from_python_subtree_" + sub_tree + ".out"


if run_phylop:
    # With unix !!! # wsl
    phylop_command = (
        "/mnt/c/Code/Github/phast/bin/phyloP --method SCORE "  
        "--subtree calJac3-saiBol1 --mode CONACC --base-by-base "
        "'/mnt/g/My Drive/Data/hg38_470Vert/phylop30/hg38.phyloP30way_named_pruned.mod' "
        "'/mnt/g/My Drive/Data/hg38_470Vert/multiz30/chr22_GL383583v2_alt.maf' "
        "> '/mnt/g/My Drive/Data/hg38_470Vert/output/phylop/phylop_from_python_subtree_calJac3-saiBol1.out'"
    )
    print("phylop command: ")
    print(phylop_command)

    # result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    subprocess.run(phylop_command, shell=True)
    print("Ran command!!!! subprocess")

    run_phylop_linux(pruned_tree_file.replace('G://', '/mnt/g/'), msa_file.replace('G://', '/mnt/g/'), output_phylop_file.replace('G://', '/mnt/g/'), sub_tree)
    print("ran using run_phylop_linux!!!")

# Run entire algorithm: determine the best tree division for acceleration and conservation for a segment
if run_phylop_recursive:
    print("Running tree partitioning!!!")
    # Example usage
#    output_file = "/mnt/g/My Drive/Data/hg38_470Vert/output/phylop/phylop_from_python_subtree_calJac3-saiBol1.out"
    phylop_data = read_phylop_output_fixed_step(output_phylop_file, plot_flag=True)

    print("read phylop output: ")
    print(type(phylop_data))
    print(len(phylop_data))
    print(phylop_data[0:10])

    print("Find best rate split!!!")
#    best_subtree, best_score = find_best_rate_split_in_tree(pruned_tree_file, msa_file, output_dir + "/phylop", method="greedy")
    best_subtree, best_score = fit_tree_rates(pruned_tree_file, msa_file, output_dir + "/phylop/chr22_GL383583v2_alt_split_tree.out", plot_tree=True)
    print("Found best rate split: ", best_subtree, best_score)
