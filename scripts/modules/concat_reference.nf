process CONCAT_REFERENCE {
    label "cpu_1"
    label "mem_1"
    label "time_1"
    
    container 'quay.io/sangerpathogens/python_graphics:1.1.5'

    input:
    path(ref)

    output:
    path(concatenated_ref)
    
    script:
    concatenated_ref = "${ref.simpleName}_concat.fasta"
    """
    concatenate_reference.py ${ref} ${concatenated_ref}
    """
}
