#!/usr/bin/env python3

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os
import sys

def main():
    if len(sys.argv) != 3:
        print("Usage: python concat_reference.py <input.fasta> <output.fasta>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    concatenated_seq = ""
    seq_records = list(SeqIO.parse(input_file, "fasta"))
    
    if not seq_records:
        print("Error: No sequences found in input file!")
        sys.exit(1)

    for record in seq_records:
        concatenated_seq += str(record.seq)

    # Extract the first fasta header (without description) for the ID
    first_record = seq_records[0]
    sequence_id = first_record.id.split()[0]

    new_record = SeqRecord(
        seq=concatenated_seq,
        id=sequence_id,
        description=f"Concatenated reference from {input_file}"
    )

    with open(output_file, "w") as out_handle:
        SeqIO.write(new_record, out_handle, "fasta")
    
    print(f"Saved concatenated sequence to {output_file} (ID: '{sequence_id}')")

if __name__ == "__main__":
    main()