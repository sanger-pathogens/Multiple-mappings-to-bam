import os
import subprocess
import sys


class Options:
    def __init__(self):
        self.ref = ''
        self.program = 'bwa'
        self.domapping = True
        self.human = False
        self.pairedend = True
        self.maxinsertsize = 1000
        self.mininsertsize = 50
        self.ssahaquality = 30
        self.maprepeats = False
        self.nomapid = 0
        self.GATK = True
        self.markdup = True
        self.detectOverlaps = False
        self.pseudosequence = True
        self.incref = True
        self.indels = True
        self.quality = 50
        self.mapq = 20
        self.depth = 8
        self.stranddepth = 3
        self.anomolous = True
        self.BAQ = True
        self.circular = True
        self.ratio = 0.8
        self.prior = 0.001
        self.call = 'c'
        self.embl = ''
        self.output = ''
        self.diroutput = ''
        self.force = False
        self.filter = '1'
        self.tabfile = False
        self.alnfile = False
        self.raxml = False
        self.model = 'GTRGAMMA'
        self.bootstrap = 100
        self.keep = False
        self.LSF = True
        self.LSFQ = 'normal'
        self.mem = 5
        self.nodes = 20
        self.dirty = False
        self.mapfiles = ''


def menu_system(option):
    run = False
    while not run:
        os.system('clear')
        print("\nINPUT OPTIONS:")
        if option.ref == '':
            print("r: Reference dna sequence:\t\tNone selected (required)")
        else:
            print("File: Directory containing the files to be mapped:\t\t" + option.mapfiles)
            print("r: Reference dna sequence:\t\t" + option.ref)
            print("p: Program:\t\t\t\t" + option.program)
            print("1: Do not remap data:\t\t\t" + str(not option.domapping))
            print("H: Mapping against human:\t\t" + str(option.human))
            print("s: Reads are single ended:\t\t" + str(not option.pairedend))
            print("i: Maximum insert size:\t\t\t" + str(option.maxinsertsize))
            print("j: Minimum insert size:\t\t\t" + str(option.mininsertsize))
            print("S: Minimum ssaha quality score:\t\t" + str(option.ssahaquality))
            print("E: Randomly map repeats:\t\t" + str(option.maprepeats))
            print("z: Minimum identity threshold:\t\t" + str(option.nomapid))
            print("G: Run GATK indel realignment:\t\t" + str(option.GATK))
            print("u: Run Mark duplicates:\t\t\t" + str(option.markdup))
            print("2: Enable read-pair overlap detection:\t" + str(option.detectOverlaps))
            print("X: Create pseudosequences:\t\t" + str(option.pseudosequence))
            print("x: Include reference in pseudosequence alignment:\t" + str(option.incref))
            print("I: Include small indels in pseudosequence alignment:\t" + str(option.indels))
            print("q: Minimum base call quality:\t\t" + str(option.quality))
            print("Q: Minimum mapping quality:\t\t" + str(option.mapq))
            print("d: Minimum number of reads matching SNP:\t" + str(option.depth))
            print("D: Minimum number of reads matching SNP per strand:\t" + str(option.stranddepth))
            print("A: Use anomolous reads in mpileup:\t" + str(option.anomolous))
            print("B: Use samtools base alignment quality option (BAQ):\t" + str(option.BAQ))
            print("c: Contigs are circular:\t\t" + str(option.circular))
            print("R: SNP/Mapping quality ratio cutoff:\t" + str(option.ratio))
            print("P: Mutation rate:\t\t\t" + str(option.prior))
            print("C: bcftools caller:\t\t\t" + option.call)
            print("e: Reference annotation:\t\t\t" + option.embl)
            print("o: Output file prefix:\t\t\t" + option.output)
            print("O: Output directory suffix:\t\t" + option.diroutput)
            print("f: Force overwrite of output files:\t" + str(option.force))
            print("F: Filter or split bam file:\t\t" + option.filter)
            print("t: Create tabfile of snps:\t\t" + str(option.tabfile))
            print("a: Create snp alignment file:\t\t" + str(option.alnfile))
            print("Y: Run phylogeny with RAxML:\t\t" + str(option.raxml))
            print("m: Model of evolution to use:\t\t" + option.model)
            print("b: Number of bootstrap replicates:\t" + str(option.bootstrap))
            print("k: If old mapping files are present, do not rerun them:\t" + str(option.keep))
            print("L: Use LSF to parallelise analyses:\t" + str(option.LSF))
            print("U: LSF queue to submit to:\t\t" + option.LSFQ)
            print("M: Amount of memory required for analysis (Gb):\t" + str(option.mem))
            print("n: Maximum number of jobs to run on nodes in parallel:\t" + str(option.nodes))
            print("y: Do not clean up temporary files:\t" + str(option.dirty))
        print("\nquit: QUIT")
        if option.ref == "":
            message = "\nPlease select an option:"
            input_list = ['r', 'Q']
        else:
            message = "\nPlease select an option or type y to run:"
            input_list = ['r', 'p', '1', 'H', 's', 'i', 'j', 'S', 'E', 'z', 'G', 'u', '2', 'X', 'x', 'I', 'q', 'Q', 'd',
                          'D', 'A', 'B', 'c', 'R', 'P', 'C', 'e', 'o', 'O', 'f', 'F', 't', 'a', 'Y', 'm', 'b', 'k', 'L',
                          'U', 'M', 'n', 'y', 'Q', 'File']
        ui = ''
        while ui not in input_list:
            ui = input(message + ' ')
        if ui == 'y':
            os.system('clear')
            run = True
        elif ui == 'File':
            option.mapfiles = input('Enter the directory containing the files to be mapped: ')
        elif ui == 'r':
            option.ref = input('Enter reference file name including path or Q to go back to the menu: ')
        elif ui == 'p':
            option.program = input('Enter program (bwa, ssaha, smalt) or Q to go back to the menu: ')
        elif ui == '1':
            option.domapping = input('Enter true or false or Q to go back to the menu: ').lower() == 'true'
        elif ui == 'H':
            option.human = input('Enter true or false or Q to go back to the menu: ').lower() == 'true'
        elif ui == 's':
            option.pairedend = input('Enter true or false or Q to go back to the menu: ').lower() == 'true'
        elif ui == 'i':
            option.maxinsertsize = int(input('Enter maximum insert size (10-10,000). Must be more than min: '))
        elif ui == 'j':
            option.mininsertsize = int(input('Enter minimum insert size (10-10,000). Must be less than max: '))
        elif ui == 'S':
            option.ssahaquality = int(input('Enter minimum ssaha quality score while mapping (ssaha only): '))
        elif ui == 'E':
            option.maprepeats = input('Enter true or false or Q to go back to the menu: ').lower() == 'true'
        elif ui == 'z':
            option.nomapid = float(input('Enter minimum identity threshold to report a mapping: '))
        elif ui == 'G':
            option.GATK = input('Enter true or false or Q to go back to the menu: ').lower() == 'true'
        elif ui == 'u':
            option.markdup = input('Enter true or false or Q to go back to the menu: ').lower() == 'true'
        elif ui == '2':
            option.detectOverlaps = input('Enter true or false or Q to go back to the menu: ').lower() == 'true'
        elif ui == 'X':
            option.pseudosequence = input('Enter true or false or Q to go back to the menu: ').lower() == 'true'
        elif ui == 'x':
            option.incref = input('Enter true or false or Q to go back to the menu: ').lower() == 'true'
        elif ui == 'I':
            option.indels = input('Enter true or false or Q to go back to the menu: ').lower() == 'true'
        elif ui == 'q':
            option.quality = int(input('Enter minimum base call quality: '))
        elif ui == 'Q':
            option.mapq = int(input('Enter minimum mapping quality: '))
        elif ui == 'd':
            option.depth = int(input('Enter minimum number of reads matching SNP: '))
        elif ui == 'D':
            option.stranddepth = int(input('Enter minimum number of reads matching SNP per strand: '))
        elif ui == 'A':
            option.anomolous = input('Enter true or false or Q to go back to the menu: ').lower() == 'true'
        elif ui == 'B':
            option.BAQ = input('Enter true or false or Q to go back to the menu: ').lower() == 'true'
        elif ui == 'c':
            option.circular = input('Enter true or false or Q to go back to the menu: ').lower() == 'true'
        elif ui == 'R':
            option.ratio = float(input('Enter SNP/Mapping quality ratio cutoff: '))
        elif ui == 'P':
            option.prior = float(input('Enter mutation rate: '))
        elif ui == 'C':
            option.call = input('Enter bcftools caller (c or m) or Q to go back to the menu: ')
        elif ui == 'e':
            option.embl = input('Enter reference annotation or Q to go back to the menu: ')
        elif ui == 'o':
            option.output = input('Enter output file prefix or Q to go back to the menu: ')
        elif ui == 'O':
            option.diroutput = input('Enter output directory suffix or Q to go back to the menu: ')
        elif ui == 'f':
            option.force = input('Enter true or false or Q to go back to the menu: ').lower() == 'true'
        elif ui == 'F':
            option.filter = input('Enter filter or split bam file (1, 2, 3, 4, 5) or Q to go back to the menu: ')
        elif ui == 't':
            option.tabfile = input('Enter true or false or Q to go back to the menu: ').lower() == 'true'
        elif ui == 'a':
            option.alnfile = input('Enter true or false or Q to go back to the menu: ').lower() == 'true'
        elif ui == 'Y':
            option.raxml = input('Enter true or false or Q to go back to the menu: ').lower() == 'true'
        elif ui == 'm':
            option.model = input(
                'Enter model of evolution to use (GTRGAMMA, GTRGAMMAI, GTRCAT, GTRMIX, GTRMIXI) or Q to go back to the menu: ')
        elif ui == 'b':
            option.bootstrap = int(input('Enter number of bootstrap replicates: '))
        elif ui == 'k':
            option.keep = input('Enter true or false or Q to go back to the menu: ').lower() == 'true'
        elif ui == 'L':
            option.LSF = input('Enter true or false or Q to go back to the menu: ').lower() == 'true'
        elif ui == 'U':
            option.LSFQ = input(
                'Enter LSF queue to submit to (normal, long, basement, hugemem) or Q to go back to the menu: ')
        elif ui == 'M':
            option.mem = int(input('Enter amount of memory required for analysis (Gb): '))
        elif ui == 'n':
            option.nodes = int(input('Enter maximum number of jobs to run on nodes in parallel: '))
        elif ui == 'y':
            option.dirty = input('Enter true or false or Q to go back to the menu: ').lower() == 'true'
        elif ui == 'quit':
            sys.exit()
    return option


def build_command(option):
    command = "nextflow run scripts/main.nf"
    if option.mapfiles is not None:
        command += f" --mapfiles {option.mapfiles}"
    if option.ref:
        command += f" --ref {option.ref}"
    if option.program:
        command += f" --program {option.program}"
    if option.domapping is not None:
        command += f" --domapping {option.domapping}"
    if option.human is not None:
        command += f" --human {option.human}"
    if option.pairedend is not None:
        command += f" --pairedend {option.pairedend}"
    if option.maxinsertsize:
        command += f" --maxinsertsize {option.maxinsertsize}"
    if option.mininsertsize:
        command += f" --mininsertsize {option.mininsertsize}"
    if option.ssahaquality:
        command += f" --ssahaquality {option.ssahaquality}"
    if option.maprepeats is not None:
        command += f" --maprepeats {option.maprepeats}"
    if option.nomapid:
        command += f" --nomapid {option.nomapid}"
    if option.GATK is not None:
        command += f" --GATK {option.GATK}"
    if option.markdup is not None:
        command += f" --markdup {option.markdup}"
    if option.detectOverlaps is not None:
        command += f" --detectOverlaps {option.detectOverlaps}"
    if option.pseudosequence is not None:
        command += f" --pseudosequence {option.pseudosequence}"
    if option.incref is not None:
        command += f" --incref {option.incref}"
    if option.indels is not None:
        command += f" --indels {option.indels}"
    if option.quality:
        command += f" --quality {option.quality}"
    if option.mapq:
        command += f" --mapq {option.mapq}"
    if option.depth:
        command += f" --depth {option.depth}"
    if option.stranddepth:
        command += f" --stranddepth {option.stranddepth}"
    if option.anomolous is not None:
        command += f" --anomolous {option.anomolous}"
    if option.BAQ is not None:
        command += f" --BAQ {option.BAQ}"
    if option.circular is not None:
        command += f" --circular {option.circular}"
    if option.ratio:
        command += f" --ratio {option.ratio}"
    if option.prior:
        command += f" --prior {option.prior}"
    if option.call:
        command += f" --call {option.call}"
    if option.embl:
        command += f" --embl {option.embl}"
    if option.output:
        command += f" --output {option.output}"
    if option.diroutput:
        command += f" --diroutput {option.diroutput}"
    if option.force is not None:
        command += f" --force {option.force}"
    if option.filter:
        command += f" --filter {option.filter}"
    if option.tabfile is not None:
        command += f" --tabfile {option.tabfile}"
    if option.alnfile is not None:
        command += f" --alnfile {option.alnfile}"
    if option.raxml is not None:
        command += f" --raxml {option.raxml}"
    if option.model:
        command += f" --model {option.model}"
    if option.bootstrap:
        command += f" --bootstrap {option.bootstrap}"
    if option.keep is not None:
        command += f" --keep {option.keep}"
    if option.LSF is not None:
        command += f" --LSF {option.LSF}"
    if option.LSFQ:
        command += f" --LSFQ {option.LSFQ}"
    if option.mem:
        command += f" --mem {option.mem}"
    if option.nodes:
        command += f" --nodes {option.nodes}"
    if option.dirty is not None:
        command += f" --dirty {option.dirty}"
    command += f" -process.echo"
    return command


def run_pipeline(option):
    command = build_command(option)
    print(command)
    subprocess.run(command, shell=True)


if __name__ == "__main__":
    options = Options()
    options = menu_system(options)
    run_pipeline(options)
