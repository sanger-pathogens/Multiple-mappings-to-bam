process HANDLE_SEQUENCES {
    container 'quay.io/ssd28/gsoc-experimental/biopython:0.0.1' // Docker container with Biopython

    input:
    path ref  // Reference FASTA file

    output:
    path("${ref}") // Updated reference file
    path("*.aln"), optional: true
    
    script:
    human = (params.human) ? "True" : "False"
    incref = (params.incref) ? "True" : "False"
    pseudosequence = (params.pseudosequence) ? "True" : "False"
    output = params.output
    """
    #!/usr/bin/env python
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    import sys
    import os

    def DoError(message):
        sys.stderr.write(message + "\\n")
        sys.exit(1)

    ref_file = "${ref}"
    incref = ${incref}
    human = ${human}
    output = "${output}"
    pseudosequence = ${pseudosequence}

    seq_records = []

    for seq_record in SeqIO.parse(open(ref_file), "fasta"):
        seq_records.append(seq_record)

    if not human:
        if len(seq_records) == 0:
            DoError("Cannot open reference fasta file!")
        else:
            SeqIO.write(seq_records, open(ref_file, "w"), "fasta")
            if incref:
                concatenated_seq = ""

                contigs = {}
                contigorder = []
                for record in seq_records:
                    contigorder.append(record.id)
                    contigs[record.id] = str(record.seq)

                keys = sorted(contigs.keys())
                for contig in contigorder:
                    concatenated_seq += contigs[contig]

                my_seq_record = SeqRecord(Seq(concatenated_seq))
                my_seq_record.id = ref_file.split("/")[-1].split(".")[0]
                my_seq_record.description = "Reference"

                # Uncomment this if you want to save a pseudosequence alignment file
                if pseudosequence:
                    SeqIO.write([my_seq_record], open(output+".aln","w"), "fasta")
            else:
                os.system("rm " + ref_file + ".aln")
    """
}
