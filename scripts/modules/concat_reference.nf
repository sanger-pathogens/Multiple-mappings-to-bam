process CONCAT_REFERENCE {
    container 'quay.io/sangerpathogens/python_graphics:1.1.5' // Docker container with Biopython

    input:
    path(ref)  // Reference FASTA file

    output:
    path(concatenated_ref) // Updated reference file
    
    script:
    concatenated_ref = "${ref.simpleName}_concat.fasta"
    """
    concatenate_reference.py ${ref} ${concatenated_ref}
    """
}
