import copy

from Bio import AlignIO, Phylo
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import tempfile
import subprocess

from utils import *
from phylo_utils import *
from subst_models import *
import os


def convert_phylip_to_maf(phylip_file, maf_file, ref_species="hg38", start_position=0):
    """
    Convert a PHYLIP alignment file to MAF format for phyloP.

    Parameters:
    - phylip_file: Path to the input PHYLIP file.
    - maf_file: Path to the output MAF file.
    - ref_species: Name of the reference species in the alignment (optional).
    - start_position: Starting position for each sequence (default is 0).

    Output:
    - MAF file written to `maf_file`.
    """
    # Read the PHYLIP alignment
    alignment = AlignIO.read(phylip_file, "phylip")

    with open(maf_file, "w") as maf_out:
        # Write MAF header
        maf_out.write("##maf version=1 scoring=none\n\n")

        # Define alignment block
        maf_out.write("a score=0.0\n")

        # Convert each sequence to MAF format
        for record in alignment:
            species = record.id
            sequence = str(record.seq)
            size = len(sequence.replace("-", ""))  # number of aligned (non-gap) bases

            # Example strand and chromosome placeholders; adjust as needed
            strand = "+"
            chromosome = f"chr1"  # Adjust to your specific chromosome if known

            # Reference start position is used for all sequences; could vary by species if known
            maf_out.write(
                f"s {species}.{chromosome:<15} {start_position:<10} {size:<10} {strand} {start_position + size:<10} {sequence}\n")

        maf_out.write("\n")  # Blank line after block



def filter_msa_by_species(msa, include_species):
    """Filter the MSA to keep only sequences for species in the tree."""
    return MultipleSeqAlignment([record for record in msa if record.id in include_species])


def filter_msa_blocks_by_species(msa_blocks, include_species):
    print("Filtering alignment by species: ")
    """Filter alignment blocks by keeping only sequences for species in the tree."""
    if isinstance(msa_blocks, str): # input as file name
        print("Load msa from file ..")
        msa_blocks = AlignIO.parse(msa_blocks, "maf")
    filtered_blocks = []
    for block in msa_blocks:
        print("Filtering block!!")
        print(block)
        filtered_records = [record for record in block if record.id in include_species]
        if filtered_records:
            filtered_block = MultipleSeqAlignment(filtered_records)
            print("Filtered: ")
            print(filtered_block)
            filtered_blocks.append(filtered_block)
    return filtered_blocks


def extract_subalignment(msa_file, start, end, output_file, include_species=None):
    """
    Extracts a sub-alignment from the given MSA file (in MAF format) between specified start and end positions.

    Parameters:
        msa_file (str): Path to the input MAF file.
        start (int): Start position for the sub-alignment.
        end (int): End position for the sub-alignment.
        output_file (str): Path to save the sub-alignment MAF file.
        include_species (list, optional): extract sequences only for these species
    """
    with open(msa_file, "r") as input_handle, open(output_file, "w") as output_handle:
        # Parse the MAF alignment file
        alignment_blocks = AlignIO.parse(input_handle, "maf")
        if include_species is not None:  # include only species present
            alignment_blocks = filter_msa_blocks_by_species(alignment_blocks, include_species)

        selected_blocks = []

        block_ctr = 0
        for block in alignment_blocks:
            block_ctr += 1


            # Check each sequence in the block
            block_in_range = False
            for record in block:
                start_pos_block = record.annotations["start"]
                end_pos_block = start_pos_block + record.annotations["size"]

                # Adjust for reverse strand if necessary
                if record.annotations["strand"] == -1:
                    start_pos_block, end_pos_block = end_pos_block, start_pos_block

                # Determine if block is within the specified range
                if min(start_pos_block, end_pos_block) <= end and max(start_pos_block, end_pos_block) >= start: #  if not (end < start_pos or start > end_pos):
#                    print("Found block in range!!, ctr=", block_ctr)
                    block_in_range = True
                    break

            # If any sequence in the block overlaps the specified range, add the block
            if block_in_range:
                updated_records = []
                new_start = max(0, start-start_pos_block)
                new_end = min(end,end_pos_block) - start_pos_block
                if new_end <= new_start:
                    continue  # skip empty blocks
                for record in block:
                    # Extract the sub-sequence using 0-based slicing
                    sub_seq = record.seq[new_start:new_end]
                    # Update the start and size annotations
                    new_start_annotation = record.annotations["start"] + new_start
                    new_size_annotation = len(sub_seq)

                    # Create a copy of the record with the updated sequence and annotations
#                    print("original record: ", record)
                    updated_record = copy.deepcopy(record)
#                    print("copy record: ", updated_record)
                    updated_record.seq = sub_seq
                    updated_record.annotations["start"] = new_start_annotation
                    updated_record.annotations["size"] = new_size_annotation
#                    print("start: ", new_start_annotation, " size=", new_size_annotation)
                    updated_records.append(updated_record)
                # Chop block!!!
                print("Found block!!! start=", start, " ; end=", end, "; start_pos=", start_pos_block, "; end_pos=", end_pos_block, " add len=", end_pos_block-start_pos_block+1,
                      " new_start=", new_start, " new_end=", new_end, " block_ctr=", block_ctr)
#                selected_blocks.append(block)
                selected_blocks.append(MultipleSeqAlignment(updated_records))

        # Write selected blocks to the output file
        AlignIO.write(selected_blocks, output_handle, "maf")
#    print(f"Sub-alignment saved to {output_file}")
    return selected_blocks


def prune_tree_by_msa_species(tree_file, msa_file):
    msa_blocks = AlignIO.parse(msa_file, "maf")
    msa_species = set()

    for block in msa_blocks:
        for record in block:
            # Extract the species name (assuming it's the prefix before a dot)
            species = record.id.split('.')[0]
            msa_species.add(species)
#            start_pos_vec.append(record.annotations["start"])
#            end_pos_vec.append(start_pos_vec[-1] + record.annotations["size"])

#    print("start: ", start_pos_vec)
#    print("end: ", end_pos_vec)

    # Extract species names in the tree
    tree = Phylo.read(tree_file, "newick")
    tree_species = {leaf.name for leaf in tree.get_terminals()}

    # Find species missing in the MSA
    missing_species_in_msa = tree_species - msa_species
    pruned_tree_file = tree_file[:-4] + "_pruned.mod"
    return extract_and_prune_model(tree_file, missing_species_in_msa, pruned_tree_file), pruned_tree_file


def filter_and_concatenate_alignment(msa_file, subtree_species, output_file=None):
    """
    Filters positions across multiple alignment blocks and concatenates them into a single block.
    Only positions with at least one non-gap nucleotide in both the subtree and its complement are retained.

    Parameters:
    - msa_file: Path to the input MSA file (in MAF format).
    - subtree_species: List of species in the subtree (matching the IDs in the alignment).
    - output_file: Path to save the concatenated filtered alignment.

    Returns:
    - good_positions: List of positions (1-based) retained in the concatenated alignment.
    """
    # Parse the MSA file as a generator of alignment blocks
    msa_blocks = AlignIO.parse(msa_file, "maf")

    # Initialize containers for the filtered sequences and good positions
    concatenated_seqs = {}
    good_positions = []
    current_position = 1
    start_positions = {}
    sequence_sizes = {}
    src_sizes = {}

#    print("filtering subtree-species=", subtree_species)
    # Process each alignment block
    for block in msa_blocks:
        alignment_length = block.get_alignment_length()

        # Split records into subtree and non-subtree species
        subtree_records = [record for record in block if record.id.split('.')[0] in subtree_species]
        other_records = [record for record in block if record.id.split('.')[0] not in subtree_species]

        # Initialize start positions for each sequence from the first block
        if not start_positions:
            for record in block:
                start_positions[record.id] = int(record.annotations.get("start", 0))
                sequence_sizes[record.id] = 0
                src_sizes[record.id] = int(record.annotations.get("srcSize", 0))

        # Filter positions in the current block
        for i in range(alignment_length):
            # Count non-gap nucleotides in the subtree and non-subtree species
            subtree_non_gaps = sum(1 for record in subtree_records if record.seq[i] != '-')
            other_non_gaps = sum(1 for record in other_records if record.seq[i] != '-')

            # print("Column: ", i, ' subtree_nt=', subtree_non_gaps, ' comp_nt=', other_non_gaps)

            # Keep the position if both subtree and non-subtree have at least one non-gap
            if max(subtree_non_gaps,other_non_gaps) > 1 and min(subtree_non_gaps, other_non_gaps) > 0:
                good_positions.append(current_position + i)

                # Append the column to each sequence in the concatenated alignment
                for record in block:
 #                   print("Update record.id: ", record.id)
                    if record.id not in concatenated_seqs:
                        concatenated_seqs[record.id] = []
                    concatenated_seqs[record.id].append(record.seq[i])
                    sequence_sizes[record.id] += 1  # Accumulate the length for this record

        current_position += alignment_length  # Update the global position for the next block

    # Create the concatenated alignment block
    concatenated_records = [SeqRecord(Seq("".join(seq)), id=seq_id, description="",
                            annotations={"start": start_positions[seq_id], "size": sequence_sizes[seq_id], "srcSize": src_sizes[seq_id]})
                            for seq_id, seq in concatenated_seqs.items()]
#    print(["start: "+ str(start_positions[seq_id]) + "size:" + str(sequence_sizes[seq_id]) for seq_id, seq in concatenated_seqs.items()])
    concatenated_alignment = MultipleSeqAlignment(concatenated_records)

    # Save the concatenated filtered alignment to the output file
    if output_file is None:
        output_file = msa_file.replace('.maf', '_filtered.maf')
    AlignIO.write(concatenated_alignment, output_file, "maf")

#    print(f"Filtered concatenated alignment saved to {output_file}")
    bad_positions = list(set(range(1, alignment_length+1)) - set(good_positions))
    return good_positions, bad_positions



def write_rate_matrix(matrix, filename):
    """
    Writes a 4x4 rate matrix to a file in the format IQ-TREE expects.

    Parameters:
    - matrix: A 4x4 rate matrix as a list of lists.
    - filename: Path to the output file.
    """
    with open(filename, 'w') as f:
        for row in matrix:
            f.write(" ".join(map(str, row)) + "\n")


def write_branch_rates(rates, filename):
    """
    Writes branch-specific rates to a file in IQ-TREE format.

    Parameters:
    - rates: A dictionary mapping branch identifiers to rates.
    - filename: Path to the output branch rates file.
    """
    with open(filename, 'w') as f:
        for branch_id, rate in rates.items():
            f.write(f"{branch_id} {rate}\n")


def parse_mode_file(mod_file):
    """
    Parses a .mod file to extract the rate matrix only.

    Parameters:
    - mod_file: Path to the model (.mod) file.

    Returns:
    - rate_matrix: 4x4 rate matrix as a list of lists.
    """
    tree = Phylo.read(mod_file, "newick")

    rate_matrix = []
    rate_matrix_section = False

    with open(mod_file, 'r') as file:
        for line in file:
#            print("read line: ", line)
            # Check if we've reached the RATE_MAT section
            if "RATE_MAT" in line:
                rate_matrix_section = True
                continue

            # Collect rate matrix values
            if rate_matrix_section:
                # End of RATE_MAT section if we reach an empty or non-numeric line
                if not line.strip() or not re.match(r'^-?\d+(\.\d+)?(\s+-?\d+(\.\d+)?)*$', line.strip()):
                    break

                # Parse a row of the rate matrix
                row = list(map(float, line.strip().split()))
                rate_matrix.append(row)


    # Ensure matrix is 4x4
    print("Rate matrix: ", rate_matrix)
    if len(rate_matrix) != 4 or any(len(row) != 4 for row in rate_matrix):
        raise ValueError("Invalid rate matrix in the .mod file. Expected a 4x4 matrix.")

    return tree, rate_matrix


def simulate_alignment_iqtree(mod_file, alignment_length, branch_rates, output_prefix="alignment"):
    """
    Simulate a multiple sequence alignment using IQ-TREE 2 based on a tree, rate matrix, and branch-specific rates.

    Parameters:
    - mod_file: Path to the Newick tree file.
    - alignment_length: Length of the simulated alignment.
    - branch_rates: A dictionary of branch-specific rates.
    - output_prefix: Prefix for the output alignment files.

    Returns:
    - output_file: Path to the generated alignment file.
    """

    def simplify_mod_file_to_iqtree(input_mod_file, output_mod_file):
        """
        Simplifies a complex .mod file to a format compatible with IQ-TREE by removing unsupported sections.

        Parameters:
        - input_mod_file: Path to the original .mod file.
        - output_mod_file: Path to save the simplified .mod file for IQ-TREE.
        """
        with open(input_mod_file, 'r') as infile, open(output_mod_file, 'w') as outfile:
            for line in infile:
                # Skip unsupported sections
                if line.startswith("ALPHABET:") or line.startswith("RATE_MAT:") or line.startswith("BACKGROUND:") or line.startswith("ORDER:") or line.startswith("SUBST_MOD:") or line.startswith("TREE:"):
                    continue

                # Write only supported model lines
                outfile.write(line)

        print(f"Simplified model file saved to {output_mod_file}")

#    simplify_mod_file_to_iqtree(mod_file, mod_file.replace('.mod', '_iqtree.mod'))

    # Step 1: Write the tree and rate matrix to a temporary file
    tree, rate_matrix = parse_mode_file(mod_file)
    # Write the tree to a temporary file
    tree_file = tempfile.NamedTemporaryFile(delete=False, mode='w', suffix=".nwk")
    Phylo.write(tree, tree_file, "newick")
    tree_file.close()

    matrix_file = tempfile.NamedTemporaryFile(delete=False, mode='w', suffix=".matrix")
    write_rate_matrix(rate_matrix, matrix_file.name)
    matrix_file.close()

    # Step 2: Write the branch rates to a temporary file
    branch_rates_file = tempfile.NamedTemporaryFile(delete=False, mode='w', suffix=".rates")
    write_branch_rates(branch_rates, branch_rates_file.name)
    branch_rates_file.close()

    rev_rates, rev_freqs = fit_gtr_model(np.array(rate_matrix))
    print("fitted reversible model rates=", rev_rates, " ; freqs=", rev_freqs)
    save_gtr_model(rev_rates, rev_freqs, mod_file.replace('.mod', '_iqtree.mod'))

    # Run IQ-TREE simulation with tree, rate matrix, and branch rates
    output_file = f"{output_prefix}.phy"
    iqtree_cmd = [
        "/mnt/c/Code/IQTree/iqtree-2.3.6-Linux-intel/bin/iqtree2",  # Update this path if necessary
        "--alisim", output_prefix,
        "--seqtype", "DNA",
        "-t", tree_file.name,  # Use the file path instead of file object
        "-m", 'GTR',  #        "-m", mod_file.replace('.mod', '_iqtree.mod'),  # Use the .mod file directly as the model input
        "--length", str(alignment_length),
#        "--ratefile", branch_rates_file.name
    ]

#    iqtree_cmd = [
#        "/mnt/c/Code/IQTree/iqtree-2.3.6-Linux-intel/bin/iqtree2",  # Update this path if necessary
#        "--alisim", output_prefix,
#        "-t", tree_file.name,  # Use the file path instead of file object
#        "-m", "GTR+G",  # Using GTR model as an example
#        "--length", str(alignment_length),
#        "--matrixfile", matrix_file.name,
#        "--ratefile", branch_rates_file.name
#    ]

    print("mod_file=", mod_file, " ; Running simulation command: ",iqtree_cmd)
    try:
        subprocess.run(iqtree_cmd, check=True)
        print(f"Simulated alignment saved to {output_file}")
    finally:
        # Clean up temporary files
        os.remove(matrix_file.name)
        os.remove(branch_rates_file.name)

    # Convert to .maf file:
#    phylip_file = "alignment.phy"  # Input PHYLIP file from IQ-TREE
    convert_phylip_to_maf(output_file, output_file.replace('.phy', '.maf'))
#    print(f"Converted MAF file saved to {maf_file}")
    # write also rates:
    write_dict_to_file(branch_rates, mod_file.replace('.mod', '_branch_rates.txt'))

    return output_file


