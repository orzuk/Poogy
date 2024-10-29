# Run phylop routines
from Bio import Phylo
import matplotlib.pyplot as plt
import os
from Bio import AlignIO

import numpy as np
from matplotlib import cm, colors
from phylo_utils import *
from phylo_plot_utils import *
from msa_utils import *
from sys import platform

# Operations to do:
reload_modules = False
plot_trees = False
prune_tree_by_msa_species = False
run_phylop = False
run_phylop_recursive = False
parse_msa = True
alignment = "multiz470" # "multiz470"  # multiz30


if reload_modules:
    import phylo_utils
    import phylo_plot_utils
    import msa_utils
    importlib.reload(phylo_utils)
    importlib.reload(msa_utils)
    importlib.reload(phylo_plot_utils)
    from phylo_utils import *
    from phylo_plot_utils import *
    from msa_utils import *


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

if alignment == "multiz470": # Run entire algorithm: determine the best tree division for acceleration and conservation for a segment
    tree_file = data_dir + "/phylop470/hg38.phyloP470way.mod"
    msa_file = data_dir + "/multiz470/PAR_CHR8_5217000_5217500.maf"  # why one block?
    phylop_str = "phylop470"
    start_pos, end_pos = 5217100, 5217800
else:
    tree_file = data_dir + "/phylop30/hg38.phyloP30way_named.mod"
    msa_file = data_dir + "/multiz30/chr22_GL383583v2_alt.maf"
    phylop_str = "phylop30"
    start_pos, end_pos = None, None  # do not extract sub-alignment

element_str = os.path.basename(msa_file)[:-4]
output_dir = data_dir + "/output"
sub_tree = "calJac3-saiBol1"
output_phylop_file = data_dir + "/output/phylop/phylop_from_python_subtree_" + sub_tree + ".out"

# Read the tree from a Newick file. Prune tree (Hard coded!)
tree = Phylo.read(tree_file, "newick")
if prune_tree_by_msa_species: # Prune tree and save
    missing_species_in_msa = ['ponAbe2', 'papAnu3', 'nasLar1', 'rhiRox1', 'rhiBie1', 'tarSyr2', 'micMur3', 'proCoq1', 'eulMac1', 'eulFla1', 'mm10', 'dasNov3']
    pruned_tree = prune_tree(tree, missing_species_in_msa)
    pruned_tree_file = tree_file[:-4] + "_pruned.mod"
    Phylo.write(pruned_tree, pruned_tree_file, "newick")
    print("Prune tree from mod file:")
    #extract_and_prune_model(tree_file, missing_species_in_msa, pruned_tree_file)
    print("READ Prune tree from mod file:")
    pruned_tree = Phylo.read(pruned_tree_file, "newick")
    # Run phylop:


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


    print("Color trees:")
    color_tree(tree, subtree_name="hg38-calJac3", output_file=data_dir + "/output/trees/output_tree_hg38-calJac3_subtree.png")
    color_tree(tree, subtree_name="calJac3-saiBol1", output_file=data_dir + "/output/trees/output_tree_calJac3-saiBol1_subtree.png")
    color_tree(pruned_tree, subtree_name="hg38-calJac3", output_file=data_dir + "/output/trees/output_pruned_tree_hg38-calJac3_subtree.png")
    color_tree(pruned_tree, subtree_name="hg38-rheMac8", output_file=data_dir + "/output/trees/output_pruned_tree_hg38-rheMac8_subtree.png")

    values_dict = dict(zip(tree.get_terminals() + tree.get_nonterminals(),
                           np.random.random(len(tree.get_terminals()) + len(tree.get_nonterminals()))))
    color_tree(tree, values=values_dict, output_file=data_dir + "/output/trees/output_tree_hg38_randcolors.png")


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


if run_phylop_recursive:
    print("Running tree partitioning!!!")
    # Example usage
#    output_file = "/mnt/g/My Drive/Data/hg38_470Vert/output/phylop/phylop_from_python_subtree_calJac3-saiBol1.out"
#    phylop_data = read_phylop_output_fixed_step(output_phylop_file, plot_flag=True)

#    print("read phylop output: ")
#    print(type(phylop_data))
#    print(len(phylop_data))
#    print(phylop_data[0:10])

    # read alignment:
    msa_blocks = AlignIO.parse(msa_file, "maf")
    msa_species = set()
    for block in msa_blocks:
        for record in block:
            # Extract the species name (assuming it's the prefix before a dot)
            species = record.id.split('.')[0]
            msa_species.add(species)

    # Extract species names in the tree
    tree_species = {leaf.name for leaf in tree.get_terminals()}

    # Find species missing in the MSA
    missing_species_in_msa = tree_species - msa_species
#    print("pruning:")
#    print(tree)
#    print(missing_species_in_msa)
#    if not tree.rooted:
#        tree.root_with_outgroup(tree.root)
#    pruned_tree = prune_and_update_tree(tree, missing_species_in_msa)
#    print("pruned!!")
    pruned_tree_file = tree_file[:-4] + "_pruned.mod"
    extract_and_prune_model(tree_file, missing_species_in_msa, pruned_tree_file)

    print("Find best rate split!!!")
    best_subtree, best_score = fit_tree_rates(pruned_tree_file, msa_file, output_dir + "/" + phylop_str + "/" + element_str + "_split_tree.out", plot_tree=True)
    print("Found best rate split: ", best_subtree, best_score)


if parse_msa:
    tree30 = Phylo.read(data_dir + "/phylop30/hg38.phyloP30way_named.mod", "newick") # TEMP HARD CODED 30
    msa_output_file = msa_file[:-4] + "_pos_" + str(start_pos) + "_" + str(end_pos) + ".maf"
    msa_selected_blocks = extract_subalignment(msa_file, start_pos, end_pos, msa_output_file)
    msa_selected_blocks2 = AlignIO.parse(msa_file, "maf")

    values_dict = dict(zip(tree30.get_terminals() + tree30.get_nonterminals(),
            np.random.random(len(tree30.get_terminals()) + len(tree30.get_nonterminals()))))
    # Now color tree with alignment !!!

    color_tree(tree30, msa=msa_selected_blocks, values=values_dict,
               output_file=data_dir + "/output/trees/output_tree_with_msa_new_filtered.png")
    color_tree(tree30, values=values_dict,
               output_file=data_dir + "/output/trees/output_tree_no_msa_new_filtered.png")
