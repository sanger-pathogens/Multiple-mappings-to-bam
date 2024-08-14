import os
import subprocess
import sys


class Options:
    # Function to initialise default values to the options
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
        self.mapfiles = []
        self.outdir = "results"

# Function for UI of CLI and input validity
def menu_system(option):
    run = False
    while not run:
        os.system('clear')
        print("\nINPUT OPTIONS:")
        if option.ref == '':
            print("r: Reference dna sequence:\t\tNone selected (required)")
        else:
            print("File: Directory containing the files to be mapped:\t\t" + str(option.mapfiles))
            print("r: Reference dna sequence:\t\t" + option.ref)
            print("outdir: Output directory name to dump all the outputs:\t\t" + options.outdir)
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
            input_list = ['r', 'QUIT']
        else:
            message = "\nPlease select an option or type y to run:"
            input_list = ['r', 'p', '1', 'H', 's', 'i', 'j', 'S', 'E', 'z', 'G', 'u', '2', 'X', 'x', 'I', 'q', 'Q', 'd',
                          'D', 'A', 'B', 'c', 'R', 'P', 'C', 'e', 'o', 'O', 'f', 'F', 't', 'a', 'Y', 'm', 'b', 'k', 'L',
                          'U', 'M', 'n', 'y', 'Q', 'File', 'outdir', 'QUIT']
        ui = ''
        while ui not in input_list:
            ui = input(message + ' ')
        if ui == 'y':
            os.system('clear')
            run = True
        elif ui == 'outdir':
            file_path = input('Enter the path of the output directory')
            options.outdir = file_path
        elif ui == 'File':
            file_path = input('Enter the path of the files to be mapped: ')
            if os.path.isfile(file_path):
                option.mapfiles.append(file_path)
            else:
                print("File does not exist")
        elif ui == 'r':
            while True:
                ref_file_path = input('Enter reference file name including path or Q to go back to the menu: ')
                print(ref_file_path)
                if os.path.isfile(ref_file_path):
                    option.ref = ref_file_path
                    break
                else:
                    print("Reference file does not exist")

        elif ui == 'p':
            while True:
                inp = input('Enter program (bwa, ssaha, smalt) or Q to go back to the menu: ').lower()
                if inp == 'q':
                    break
                elif inp in ['bwa', 'ssaha', 'smalt']:
                    option.program = inp
                    break
                else :
                    print("Invalid program")

        elif ui == '1':
            while True:
                inp = input('Enter true or false or Q to go back to the menu: ').lower()
                if inp == 'q':
                    break
                elif inp in ['true', 'false']:
                    option.domapping = inp == 'true'
                    break
                else :
                    print("Invalid input")
        elif ui == 'H':
            while True:
                inp = input('Enter true or false or Q to go back to the menu: ').lower()
                print(inp)
                if inp == 'q':
                    break
                elif inp in ['true', 'false']:
                    option.human = inp == 'true'
                    break
                else :
                    print("Invalid input")
        elif ui == 's':
            while True:
                inp = input('Enter true or false or Q to go back to the menu: ').lower()
                print(inp)
                if inp == 'q':
                    break
                elif inp in ['true', 'false']:
                    option.pairedend = inp == 'true'
                    break
                else :
                    print("Invalid input")
        elif ui == 'i':
            while True:
                inp = input('Enter maximum insert size (10-10,000) or Q to go back to the menu. Must be more than min: ').lower()
                if inp == 'q':
                    break
                elif int(inp) >= 10 and int(inp) <= 10000 and int(inp) > int(option.mininsertsize):
                    option.maxinsertsize = inp
                    break
                else :
                    print("Invalid input")
            
        elif ui == 'j':
            while True:
                inp = input('Enter minimum insert size (10-10,000) or Q to go back to the menu. Must be less than max: ').lower()
                if inp == 'q':
                    break
                elif int(inp) >= 10 and int(inp) <= 10000 and int(inp) < int(option.maxinsertsize):
                    option.mininsertsize = inp
                    break
                else :
                    print("Invalid input")
        elif ui == 'S':
            if option.program != 'ssaha':
                print("ssaha quality score is only applicable for ssaha program. Press any key to continue.")
                inp = input()
            while True:
                inp = input('Enter minimum ssaha quality score while mapping (ssaha only) or Q to go back to the menu: ').lower()
                if inp == 'q':
                    break
                elif int(inp) >= 0:
                    option.ssahaquality = inp
                    break
                else :
                    print("Invalid input")
        elif ui == 'E':
            while True:
                inp = input('Enter true or false or Q to go back to the menu: ').lower()
                print(inp)
                if inp == 'q':
                    break
                elif inp in ['true', 'false']:
                    option.maprepeats = inp == 'true'
                    break
                else :
                    print("Invalid input")
        elif ui == 'z':
            option.nomapid = float(input('Enter minimum identity threshold to report a mapping: '))
        elif ui == 'G':
            while True:
                inp = input('Enter true or false or Q to go back to the menu: ').lower()
                print(inp)
                if inp == 'q':
                    break
                elif inp in ['true', 'false']:
                    option.GATK = inp == 'true'
                    break
                else :
                    print("Invalid input")
        elif ui == 'u':
            while True:
                inp = input('Enter true or false or Q to go back to the menu: ').lower()
                print(inp)
                if inp == 'q':
                    break
                elif inp in ['true', 'false']:
                    option.markdup = inp == 'true'
                    break
                else :
                    print("Invalid input")
        elif ui == '2':
            while True:
                inp = input('Enter true or false or Q to go back to the menu: ').lower()
                print(inp)
                if inp == 'q':
                    break
                elif inp in ['true', 'false']:
                    option.detectOverlaps = inp == 'true'
                    break
                else :
                    print("Invalid input")
        elif ui == 'X':
            while True:
                inp = input('Enter true or false or Q to go back to the menu: ').lower()
                print(inp)
                if inp == 'q':
                    break
                elif inp in ['true', 'false']:
                    option.pseudosequence = inp == 'true'
                    break
                else :
                    print("Invalid input")
        elif ui == 'x':
            while True:
                inp = input('Enter true or false or Q to go back to the menu: ').lower()
                print(inp)
                if inp == 'q':
                    break
                elif inp in ['true', 'false']:
                    option.incref = inp == 'true'
                    break
                else :
                    print("Invalid input")
        elif ui == 'I':
            while True:
                inp = input('Enter true or false or Q to go back to the menu: ').lower()
                print(inp)
                if inp == 'q':
                    break
                elif inp in ['true', 'false']:
                    option.indels = inp == 'true'
                    break
                else :
                    print("Invalid input")
        elif ui == 'q':
            while True:
                inp = input('Enter minimum base call quality: ').lower()
                if inp == 'q':
                    break
                elif int(inp) > 0 and int(inp) < 100:
                    option.quality = inp
                    break
                else :
                    print("The value should be in range [1, 99]")
        elif ui == 'Q':
            while True:
                inp = input('Enter minimum mapping quality: ').lower()
                if inp == 'q':
                    break
                elif int(inp) > 0 and (int(inp) < 30 or option.program != 'bwa') and (int(inp) < 60 or option.program == 'bwa'):
                    option.mapq = inp
                    break
                else :
                    if option.program == 'bwa':
                        print("Mapping quality score (-Q) must be between 0 and 30 for BWA.")
                    else:
                        print(f"Mapping quality score (-Q) must be between 0 and 60 for {option.program}")
        elif ui == 'd':
            while True:
                inp = input('Enter minimum number of reads matching SNP: ').lower()
                if inp == 'q':
                    break
                elif int(inp) >= 0:
                    option.depth = inp
                    break
                else :
                    print("Invalid input")
        elif ui == 'D':
            while True:
                inp = input('Enter minimum number of reads matching SNP per strand: ').lower()
                if inp == 'q':
                    break
                elif int(inp) >= 0:
                    option.stranddepth = inp
                    break
                else :
                    print("Invalid input")
        elif ui == 'A':
            while True:
                inp = input('Enter true or false or Q to go back to the menu: ').lower()
                print(inp)
                if inp == 'q':
                    break
                elif inp in ['true', 'false']:
                    option.anomolous = inp == 'true'
                    break
                else :
                    print("Invalid input")
        elif ui == 'B':
            while True:
                inp = input('Enter true or false or Q to go back to the menu: ').lower()
                print(inp)
                if inp == 'q':
                    break
                elif inp in ['true', 'false']:
                    option.BAQ = inp == 'true'
                    break
                else :
                    print("Invalid input")
        elif ui == 'c':
            while True:
                inp = input('Enter true or false or Q to go back to the menu: ').lower()
                print(inp)
                if inp == 'q':
                    break
                elif inp in ['true', 'false']:
                    option.circular = inp == 'true'
                    break
                else :
                    print("Invalid input")
        elif ui == 'R':
            while True:
                inp = float (input('Enter SNP/Mapping quality ratio cutoff: '))
                if inp>1 or inp<0:
                    print("SNP/site mapping quality ratio cutoff (-R) must be between 0 and 1")
                else:
                    option.ratio = inp
        elif ui == 'P':
            while True:
                inp = float (input('Enter mutation rate: '))
                if inp>1 or inp<0:
                    print("Estimated mutation rate (-P) must be between 0 and 1")
                else:
                    option.prior = inp
        elif ui == 'C':
            while True:
                inp = input('Enter bcftools caller (c or m) or Q to go back to the menu: ').lower()
                if inp == 'q':
                    break
                elif inp in ['c', 'm']:
                    option.call = inp
                    break
                else :
                    print("Invalid input")
        elif ui == 'e':
            inp = input('Enter reference annotation or Q to go back to the menu: ')
            if inp != 'Q':
                option.embl = inp
        elif ui == 'o':
            inp = input('Enter output file prefix or Q to go back to the menu: ')
            if inp != 'Q':
                option.output = inp
        elif ui == 'O':
            inp = input('Enter output directory suffix or Q to go back to the menu: ')
            if inp != 'Q':
                option.diroutput = inp
        elif ui == 'f':
            while True:
                inp = input('Enter true or false or Q to go back to the menu: ').lower()
                if inp == 'q':
                    break
                elif inp in ['true', 'false']:
                    option.force = inp == 'true'
                    break
                else :
                    print("Invalid input")
        elif ui == 'F':
            while True:
                inp = input('Enter filter or split bam file (1, 2, 3, 4, 5) or Q to go back to the menu: ').lower()
                print(inp)
                if inp == 'q':
                    break
                elif inp in ['1', '2', '3', '4', '5']:
                    option.filter = inp
                    break
                else :
                    print("Invalid input")
        elif ui == 't':
            while True:
                inp = input('Enter true or false or Q to go back to the menu: ').lower()
                if inp == 'q':
                    break
                elif inp in ['true', 'false']:
                    option.tabfile = inp == 'true'
                    break
                else :
                    print("Invalid input")
        elif ui == 'a':
            while True:
                inp = input('Enter true or false or Q to go back to the menu: ').lower()
                if inp == 'q':
                    break
                elif inp in ['true', 'false']:
                    option.alnfile = inp == 'true'
                    break
                else :
                    print("Invalid input")
        elif ui == 'Y':
            while True:
                inp = input('Enter true or false or Q to go back to the menu: ').lower()
                if inp == 'q':
                    break
                elif inp in ['true', 'false']:
                    option.raxml = inp == 'true'
                    break
                else :
                    print("Invalid input")
        elif ui == 'm':
            while True:
                inp = input('Enter true or false or Q to go back to the menu: ')
                if inp == 'Q':
                    break
                elif inp in ['GTRGAMMA', 'GTRGAMMAI', 'GTRCAT', 'GTRMIX', 'GTRMIXI']:
                    option.force = inp
                    break
                else :
                    print("Invalid input")
        elif ui == 'b':
            inp = input('Enter number of bootstrap replicates: ')
            if inp == 'Q':
                continue
            elif int(inp) >= 0:
                option.bootstrap = inp
                break
        elif ui == 'k':
            while True:
                inp = input('Enter true or false or Q to go back to the menu: ').lower()
                if inp == 'q':
                    break
                elif inp in ['true', 'false']:
                    option.keep = inp == 'true'
                    break
                else :
                    print("Invalid input")
        elif ui == 'L':
            while True:
                inp = input('Enter true or false or Q to go back to the menu: ').lower()
                if inp == 'q':
                    break
                elif inp in ['true', 'false']:
                    option.LSF = inp == 'true'
                    break
                else :
                    print("Invalid input")
        elif ui == 'U':
            while True:
                inp = input('Enter LSF queue to submit to (normal, long, basement, hugemem) or Q to go back to the menu: ').lower()
                if inp == 'q':
                    break
                elif inp in ['normal', 'long', 'basement', 'hugemem']:
                    option.LSFQ = inp
                    break
                else :
                    print("Invalid input")
        elif ui == 'M':
            inp = input('Enter amount of memory required for analysis (Gb): ')
            if inp == 'Q':
                continue
            if inp >= 0 and inp <= 30:
                option.mem = inp
            else:
                print("Memory requirement (-M) must be between 0 and 30Gb")

        elif ui == 'n':
            inp = input('Enter maximum number of jobs to run on nodes in parallel: ')
            if inp == 'Q':
                continue
            if inp >=0:
                option.nodes = inp
        elif ui == 'y':
            while True:
                inp = input('Enter true or false or Q to go back to the menu: ').lower()
                if inp == 'q':
                    break
                elif inp in ['true', 'false']:
                    option.dirty = inp == 'true'
                    break
                else :
                    print("Invalid input")
        elif ui == 'QUIT':
            sys.exit()
    return option


def build_command(option):
    command = "nextflow run scripts/main.nf"
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
    if len(option.mapfiles) > 0:
        command += f" --mapfiles \""
        for file in option.mapfiles:
            command += f"{file},"
        command = command[:-1]
        command += "\""
    command += f" -process.echo"
    return command



def DoError (errorString, options):
    print ("\nError: ", errorString)
    print ("Press any key to continue.")
    inp = input()
    options = menu_system (options)
    return options



def validate_inputs(options):

    if options.output == '':
        options.output = options.ref.split("/")[-1].split(".")[0]+"_"+options.program

    while options.force==False and os.path.isfile(options.output+".aln"):
        output = ""
        output = input ("\nOutput files with chosen prefix already exist.\n\nWould you like to overwrite (o), choose a new output file prefix (n) or quit (Q): ")
        if output == 'Q':
            sys.exit()
        elif output == "o":
            options.force=True
            break
        elif output == "n":
            options.output = input("Enter a new output file prefix: ")
    if options.ref == '':
        options = DoError("No reference DNA file (-r) selected!", options)
    if not options.mapfiles:
        options = DoError("No input files selected!", options)
    options.program=options.program.upper()
    if options.mem == 0:
        options.mem = 2
    return options
    

def run_pipeline(option):
    command = build_command(option)
    print(command)
    subprocess.run(command, shell=True)


if __name__ == "__main__":
    options = Options()
    options = menu_system(options)
    options = validate_inputs(options)
    run_pipeline(options)
