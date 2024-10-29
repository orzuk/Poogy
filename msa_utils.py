from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment


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

#                print("start=", start, " ; end=", end, "; start_pos=", start_pos_block, "; end_pos=", end_pos_block)
                # Determine if block is within the specified range
                if min(start_pos_block, end_pos_block) <= end and max(start_pos_block, end_pos_block) >= start: #  if not (end < start_pos or start > end_pos):
#                    print("Found block in range!!, ctr=", block_ctr)
                    block_in_range = True
                    break

            # If any sequence in the block overlaps the specified range, add the block
            if block_in_range:
                selected_blocks.append(block)

        # Write selected blocks to the output file
        AlignIO.write(selected_blocks, output_handle, "maf")
    print(f"Sub-alignment saved to {output_file}")
    return selected_blocks


