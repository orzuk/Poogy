# Run phylop routines
from Bio import Phylo
import matplotlib.pyplot as plt
import os
from Bio import AlignIO
import importlib
import subprocess
import numpy as np
from matplotlib import cm, colors
from phylo_utils import *
from phylo_plot_utils import *
from msa_utils import *
from sys import platform
import argparse

# Create an argument parser object
parser = argparse.ArgumentParser(description="Run tree rate splitting with optional arguments.")

# Add arguments with their default values
parser.add_argument('--alignment', type=str, default='multiz470', help='Path to the tree file')
parser.add_argument('--k', type=int, default=10, help='Kmer length')
parser.add_argument('--run_mode', type=str, default='plot_trees', help='Mode to run the script')  # Here decide what to do
parser.add_argument('--msa_file', type=str, default='', help='Alignment File')  # Here decide what to do
parser.add_argument('--tree_file', type=str, default='', help='Model file with phylogenetic tree')  # Here decide what to do
parser.add_argument('--branch_rates_file', type=str, default='', help='File with branch rates')
parser.add_argument('--pos', type=str, default='', help='Genomic coordinates')  # Here decide what to do
parser.add_argument('--verbose', action='store_true', help='Enable verbose mode')  # Example for boolean operation
parser.add_argument('--h', dest='h', action='store_true', help='help')  # Example for boolean operation
parser.add_argument('--no-h', dest='h', action='store_false', help="A boolean flag that is False if provided.")



# Parse the command-line arguments
args = parser.parse_args()

# Operations to do:
reload_modules = True
plot_trees = args.run_mode == "plot_trees"
prune_tree_by_msa = args.run_mode == "prune_tree_by_msa"
run_phylop = args.run_mode == "run_phylop"
run_phylop_best_split = args.run_mode == "run_phylop_best_split"
run_phylop_recursive = args.run_mode == "plot_phylop_recursive"
parse_msa = args.run_mode == "parse_msa"
simulate_msa = args.run_mode == "simulate_msa"
help_mode = (args.run_mode == "help") or args.h
alignment = args.alignment  # "multiz470" # "multiz470"  # multiz30
try:
    alignment_length = args.k  # Length of the alignment
except ValueError:
    alignment_length = 100  # Length of the alignment

if help_mode:
    print("Usage: python3 run_phylop --alignemt=<alignment-file>, --k=<kmer-length> --run_mode=<operation-to-run>")
    sys.exit()

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

if alignment == "multiz30":  # Run entire algorithm: determine the best tree division for acceleration and conservation for a segment
    tree_file = data_dir + "/phylop470/hg38.phyloP470way.mod"
    msa_file = data_dir + "/multiz470/PAR_CHR8_5217000_5217500.maf"  # why one block?
    phylop_str = "phylop470"
    start_pos, end_pos = 5217100, 5217800
else:
    tree_file = data_dir + "/phylop30/hg38.phyloP30way_named.mod"
    msa_file = data_dir + "/multiz30/chr22_GL383583v2_alt.maf"
    phylop_str = "phylop30"
    start_pos, end_pos = 47, 82  # do not extract sub-alignment  # 12-mer
if args.pos != "":  # extract positions from input (override defaults)
    start_pos, end_pos = [int(s) for s in args.pos.split("-")]

msa_sim_output_file = "sim/simulated_alignment.maf"
if args.msa_file != "":    # override alignment stuff
    msa_file = args.msa_file
    msa_sim_output_file = msa_file
    start_pos, end_pos = None, None  # don't allow start/end pos for now
if args.tree_file != "":
    tree_file = args.tree_file
    start_pos, end_pos = None, None


if start_pos is None:
    pos_str = ""
else:
    pos_str = "_" + str(start_pos) + "_" + str(end_pos)
element_str = os.path.basename(msa_file)[:-4] + pos_str
output_dir = data_dir + "/output"
sub_tree = "calJac3-saiBol1"
output_phylop_file = data_dir + "/output/phylop/phylop_from_python_subtree_" + sub_tree + ".out"

# Read the tree from a Newick file. Prune tree (Hard coded!)
tree = Phylo.read(tree_file, "newick")
if prune_tree_by_msa: # Prune tree and save
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
    values_dict2 = dict(zip(tree.get_terminals() + tree.get_nonterminals(),
                           np.random.random(len(tree.get_terminals()) + len(tree.get_nonterminals()))))

    color_tree(tree, values=values_dict, output_file=data_dir + "/output/trees/output_tree_hg38_1st_randcolors.png")
    color_tree(tree, values=values_dict2, output_file=data_dir + "/output/trees/output_tree_hg38_2nd_randcolors.png")
    color_tree(tree, values=[values_dict, values_dict2], output_file=data_dir + "/output/trees/output_tree_hg38_BOTH_randcolors.png")


if run_phylop:  # just run on a simple split of the tree into two
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


if run_phylop_best_split:  # find automatically the best split into two complementary sub-trees
    print("Running tree partitioning to two sub-trees of different rates:")
    print(tree_file, msa_file)
    pruned_tree, pruned_tree_file = prune_tree_by_msa_species(tree_file, msa_file)  # get msa species in tree
    if start_pos is not None:
        use_msa_file = msa_file.replace('.maf', '_start_' + str(start_pos) + "_" + str(end_pos) + ".maf" )
        extract_subalignment(msa_file, start_pos, end_pos, use_msa_file)
    else:  # no need for extracting
        use_msa_file = msa_file
    output_file = output_dir + "/" + phylop_str + "/" + element_str + "_split_tree.out"
    if args.branch_rates_file != "":  # read true branch rates
        true_branch_rates = read_dict_from_file(args.branch_rates_file)
    else:
        true_branch_rates = None
    best_subtree, best_score = fit_two_subtree_rates(pruned_tree_file, use_msa_file, output_file, true_branch_rates=true_branch_rates,
                                                     plot_tree=True, method="brute-force", phylop_score="LRT")
    print("Found best rate split brute-force: ", best_subtree.root.name, best_score)
    best_subtree_binary, best_score_binary = fit_two_subtree_rates(pruned_tree_file, use_msa_file, output_file, true_branch_rates=true_branch_rates,
                                                                   plot_tree=True, method="binary", phylop_score="LRT")
    print("Found best rate split binary: ", best_subtree_binary.root.name, best_score_binary)


# split recursively into multiple subtrees with different rates
if run_phylop_recursive:
    print("Run recursive!!! ")
    # T.B.D.


if parse_msa:  # read and extract sub-alignment
    tree30 = Phylo.read(data_dir + "/phylop30/hg38.phyloP30way_named.mod", "newick") # TEMP HARD CODED 30WAY !!
#    msa_file = file_name_to_unix(msa_file)
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


if simulate_msa:
    # Example usage
#    mod_file = "example.mod"  # Path to the model file
    branch_rates = subtree_to_branch_rates(tree, 'hg38-rheMac8', subtree_rate=2.0, default_rate=1.0)  # acceleration
    simulate_alignment_iqtree(tree_file, alignment_length, branch_rates, output_prefix=msa_sim_output_file[:-4]) # remove suffix

    # Example usage
#    tree_file = "example_tree.nwk"  # Path to the Newick tree file
#    alignment_length = 1000  # Length of the alignment
#    rates = [0.1, 0.2, 0.3, 0.1, 0.5]  # Example rates for each branch
#    rate_matrix = np.array([
#        [0.0, 0.2, 0.3, 0.5],
#        [0.1, 0.0, 0.4, 0.5],
#        [0.2, 0.3, 0.0, 0.5],
#        [0.3, 0.3, 0.4, 0.0]
#    ])  # Example 4x4 rate matrix for GTR


################################################################
########### OLD STUFF/JUNK BELOW ###############################
################################################################
    # Getting back the objects:
#with open('tree_to_plot.pkl', 'rb') as f:  # Python 3: open(..., 'rb')
#        tree, subtree_rates, output_plot_tree_file = pickle.load(f)
#with open('tree_with_sub.pkl', 'rb') as f:  # Python 3: open(..., 'rb')
#    tree, best_subtree = pickle.load(f)

# with open('bad_msa_phylop_scores.pkl', 'wb') as f:  # Python 3: open(..., 'wb')
# with open('bad_msa_phylop_scores.pkl', 'rb') as f:  # Python 3: open(..., 'rb')
#    tree_file, msa_file, output_plot_tree_file, best_subtree = pickle.load(f)
################################################################
########### OLD STUFF/JUNK ABPVE ###############################
################################################################


#msa_file = "G:\My Drive\Data\hg38_470Vert\multiz30\chr22_GL383583v2_alt_start_47_82.maf"
#debug_msa_blocks = AlignIO.parse(msa_file, "maf")