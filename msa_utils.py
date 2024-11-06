import copy

from Bio import AlignIO, Phylo
from Bio.Align import MultipleSeqAlignment
from phylo_utils import *


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
                    print("original record: ", record)
                    updated_record = copy.deepcopy(record)
                    print("copy record: ", updated_record)
                    updated_record.seq = sub_seq
                    print("copy2 record: ", updated_record)
                    updated_record.annotations["start"] = new_start_annotation
                    print("copy3 record: ", updated_record)
                    updated_record.annotations["size"] = new_size_annotation
                    print("copy4 record: ", updated_record)

                    print("add update record: ", updated_record)
                    print("start: ", new_start_annotation, " size=", new_size_annotation)
                    updated_records.append(updated_record)
                # Chop block!!!
                print("Found block!!! start=", start, " ; end=", end, "; start_pos=", start_pos_block, "; end_pos=", end_pos_block, " add len=", end_pos_block-start_pos_block+1,
                      " new_start=", new_start, " new_end=", new_end, " block_ctr=", block_ctr)
#                selected_blocks.append(block)
                selected_blocks.append(MultipleSeqAlignment(updated_records))

        # Write selected blocks to the output file
        AlignIO.write(selected_blocks, output_handle, "maf")
    print(f"Sub-alignment saved to {output_file}")
    return selected_blocks


def prune_tree_by_msa_species(tree_file, msa_file):
    print("Start prunning!")
    msa_blocks = AlignIO.parse(msa_file, "maf")
    msa_species = set()

    print("Run over blocks!")
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
    print("Read tree!")
    tree = Phylo.read(tree_file, "newick")
    tree_species = {leaf.name for leaf in tree.get_terminals()}

    # Find species missing in the MSA
    missing_species_in_msa = tree_species - msa_species
    pruned_tree_file = tree_file[:-4] + "_pruned.mod"
    print("Extract and prune model!")
    return extract_and_prune_model(tree_file, missing_species_in_msa, pruned_tree_file), pruned_tree_file


