import sys
import random
from math import log, exp
import csv
from itertools import product

def create_aa_dict():
    aa_dict = {"L": "leu",
               "P": "pro",
               "H": "his",
               "Q": "gln",
               "R": "arg",
               "I": "ile",
               "M": "met",
               "T": "thr",
               "N": "asn",
               "K": "lys",
               "S": "ser",
               "V": "val",
               "A": "ala",
               "D": "asp",
               "E": "glu",
               "G": "gly",
               "F": "phe",
               "Y": "tyr",
               "C": "cys",
               "W": "trp",
               "X": "stp"}

    return aa_dict


def create_aa_counter():
    aa_counter = {"L": 0,
                  "P": 0,
                  "H": 0,
                  "Q": 0,
                  "R": 0,
                  "I": 0,
                  "M": 0,
                  "T": 0,
                  "N": 0,
                  "K": 0,
                  "S": 0,
                  "V": 0,
                  "A": 0,
                  "D": 0,
                  "E": 0,
                  "G": 0,
                  "F": 0,
                  "Y": 0,
                  "C": 0,
                  "W": 0,
                  "X": 0}
    return aa_counter


def create_codon_dict():
    codon_dict = {"AAA": "K",
                  "AAC": "N",
                  "AAG": "K",
                  "AAT": "N",
                  "ACA": "T",
                  "ACC": "T",
                  "ACG": "T",
                  "ACT": "T",
                  "AGA": "R",
                  "AGC": "S",
                  "AGG": "R",
                  "AGT": "S",
                  "ATA": "I",
                  "ATC": "I",
                  "ATG": "M",
                  "ATT": "I",
                  "CAA": "Q",
                  "CAC": "H",
                  "CAG": "Q",
                  "CAT": "H",
                  "CCA": "P",
                  "CCC": "P",
                  "CCG": "P",
                  "CCT": "P",
                  "CGA": "R",
                  "CGC": "R",
                  "CGG": "R",
                  "CGT": "R",
                  "CTA": "L",
                  "CTC": "L",
                  "CTG": "L",
                  "CTT": "L",
                  "GAA": "E",
                  "GAC": "D",
                  "GAG": "E",
                  "GAT": "D",
                  "GCA": "A",
                  "GCC": "A",
                  "GCG": "A",
                  "GCT": "A",
                  "GGA": "G",
                  "GGC": "G",
                  "GGG": "G",
                  "GGT": "G",
                  "GTA": "V",
                  "GTC": "V",
                  "GTG": "V",
                  "GTT": "V",
                  "TAA": "X",
                  "TAC": "Y",
                  "TAG": "X",
                  "TAT": "Y",
                  "TCA": "S",
                  "TCC": "S",
                  "TCG": "S",
                  "TCT": "S",
                  "TGA": "X",
                  "TGC": "C",
                  "TGG": "W",
                  "TGT": "C",
                  "TTA": "L",
                  "TTC": "F",
                  "TTG": "L",
                  "TTT": "F"}

    return codon_dict


def create_codon_counter():
    codon_counter = {"AAA": 0,
                     "AAC": 0,
                     "AAG": 0,
                     "AAT": 0,
                     "ACA": 0,
                     "ACC": 0,
                     "ACG": 0,
                     "ACT": 0,
                     "AGA": 0,
                     "AGC": 0,
                     "AGG": 0,
                     "AGT": 0,
                     "ATA": 0,
                     "ATC": 0,
                     "ATG": 0,
                     "ATT": 0,
                     "CAA": 0,
                     "CAC": 0,
                     "CAG": 0,
                     "CAT": 0,
                     "CCA": 0,
                     "CCC": 0,
                     "CCG": 0,
                     "CCT": 0,
                     "CGA": 0,
                     "CGC": 0,
                     "CGG": 0,
                     "CGT": 0,
                     "CTA": 0,
                     "CTC": 0,
                     "CTG": 0,
                     "CTT": 0,
                     "GAA": 0,
                     "GAC": 0,
                     "GAG": 0,
                     "GAT": 0,
                     "GCA": 0,
                     "GCC": 0,
                     "GCG": 0,
                     "GCT": 0,
                     "GGA": 0,
                     "GGC": 0,
                     "GGG": 0,
                     "GGT": 0,
                     "GTA": 0,
                     "GTC": 0,
                     "GTG": 0,
                     "GTT": 0,
                     "TAA": 0,
                     "TAC": 0,
                     "TAG": 0,
                     "TAT": 0,
                     "TCA": 0,
                     "TCC": 0,
                     "TCG": 0,
                     "TCT": 0,
                     "TGA": 0,
                     "TGC": 0,
                     "TGG": 0,
                     "TGT": 0,
                     "TTA": 0,
                     "TTC": 0,
                     "TTG": 0,
                     "TTT": 0}

    return codon_counter


def create_aa_syn_dict():
    aa_syn_counter = create_aa_counter()
    codon_dict = create_codon_dict()

    for aa in aa_syn_counter:
        for codon in codon_dict:
            if aa == codon_dict[codon]:
                aa_syn_counter[aa] += 1

    return aa_syn_counter


def check_seq(seq, name=""):
    seq = seq.upper()

    base_dict = ["A", "C", "G", "T", "U", "N", "-", "M", "R", "W", "S", "Y", "K", "V", "H", "D", "B"]

    for base in seq:
        if base not in base_dict:
            print("Error - check_seq - abnormal character in sequence [converted to N] = " + base + " seq = " + name)
            seq = seq.replace(base, "N")

    return seq


def remove_gaps(seq):

    return seq.upper().replace("-", "")


def replace_ambiguities(seq):
    seq = seq.upper()
    base_dict = ["A", "C", "G", "T", "U", "N", "-"]
    no_amb_seq = ''

    for base in seq:
        if base in base_dict:
            no_amb_seq += base
        else:
            no_amb_seq += 'N'

    return no_amb_seq


def correct_seq(seq):
    seq = replace_ambiguities(seq)
    seq = rna_to_dna(seq)
    seq = remove_gaps(seq)

    return seq


def count_bases(seq):
    seq = seq.upper()

    base_dict = {
        "A": 0, "C": 0, "G": 0, "T": 0, "N": 0
    }

    for base in seq:
        if base in base_dict:
            base_dict[base] += 1
        else:
            base_dict["N"] += 1

    return base_dict


def count_all_bases(seq):
    seq = seq.upper()

    base_dict = {
        "A": 0, "C": 0, "G": 0, "T": 0, "U": 0, "N": 0, "-": 0,
        "M": 0, "R": 0, "W": 0, "S": 0, "Y": 0, "K": 0, "V": 0, "H": 0, "D": 0, "B": 0,
        "Errors": 0
    }

    for base in seq:
        if base in base_dict:
            base_dict[base] += 1
        else:
            base_dict["Errors"] += 1

    return base_dict


def reverse_complement(seq, name=""):
    seq = complement(reverse(seq.upper()), name)

    return seq


def complement(seq, name=""):
    seq = seq.upper()
    c_seq = ''

    comp_dict = {
        "A": "T", "C": "G", "G": "C", "T": "A", "U": "A", "N": "N", "-": "-",
        "M": "K", "R": "Y", "W": "W", "S": "S", "Y": "R", "K": "M", "V": "B", "H": "D", "D": "H", "B": "V"
    }

    for base in seq:
        if base in comp_dict:
            c_seq += comp_dict[base]
        else:
            c_seq += 'N'
            print("Error - complement - abnormal character in sequence [converted to N] = " + base + " seq = " + name)

    return c_seq


def reverse(seq):

    return seq[::-1]


def rna_to_dna(seq):

    return seq.upper().replace("U", "T")


def dna_to_rna(seq):

    return seq.upper().replace("T", "U")


def count_dinucleotides(seq):
    seq = correct_seq(seq)

    dinuc_dict = {
        "AA": 0, "AC": 0, "AG": 0, "AT": 0,
        "CA": 0, "CC": 0, "CG": 0, "CT": 0,
        "GA": 0, "GC": 0, "GG": 0, "GT": 0,
        "TA": 0, "TC": 0, "TG": 0, "TT": 0,
        "NN": 0
    }

    for dinuc in dinuc_dict:
        dinuc_dict[dinuc] = seq.count(dinuc)

    return dinuc_dict


def count_lines_slower(filename):
    line_count = 0
    with open(filename) as file_handler:
        for line in file_handler:
            line_count += 1

    return line_count


def count_lines(filename):
    with open(filename) as f:
        for i, l in enumerate(f):
            pass

    return i + 1


def count_fasta(filename):
    seq_count = 0

    with open(filename) as file_handler:
        for line in file_handler:
            if line.find(">") == 0:
                seq_count += 1
            if line.rfind(">") > 0:
                print("Error - count_fasta - '>' symbol found mid sequence [but I'm not exiting]: " + line)

    return seq_count


def count_fastq(filename):
    seq_count = 0
    fastq_count = 0

    with open(filename) as file_handler:
        for line in file_handler:
            line = line.rstrip()
            fastq_count += 1

            if fastq_count == 1:
                read = []

            read.append(line)

            if fastq_count == 4:
                check_read(read)
                fastq_count = 0
                seq_count += 1

    if fastq_count != 0:
        print("Error - count_fastq - incomplete FASTQ record in file (" + str(seq_count) + "/" + str(fastq_count) + "): " + filename)

    return seq_count


def import_fastq(filename):
    seq_count = 0
    fastq_count = 0
    seq_reads = []

    with open(filename) as file_handler:
        for line in file_handler:
            line = line.rstrip()
            fastq_count += 1

            if fastq_count == 1:
                read = []

            read.append(line)

            if fastq_count == 4:
                check_read(read)
                seq_reads.append(get_readname_stub(read[0]) + read)
                fastq_count = 0
                seq_count += 1

    print("Reads imported from " + filename + " = " + str(seq_count))

    return seq_reads


def check_read(read):
    read_correct = True

    if read[0][0] != '@':
        read_correct = False
        print("Error - check_read - first line in FASTQ record does not begin with @: " + read[0])
    if read[2][0] != '+':
        print("Error - check_read - third line in FASTQ record does not begin with +: " + read[0] + " " + read[2])
        read_correct = False
    if len(read[1]) != len(read[3]):
        print("Error - check_read - seqLength(" + str(read[1]) + ") != qualLength(" + str(read[3]) + "): " + read[0] + " : " + read[1] + " : " + read[3])
        read_correct = False

    return read_correct


def check_pair(read1, read2):
    read1_correct = check_read(read1[1:])
    read2_correct = check_read(read2[1:])

    if read1[0] != read2[0]:
        read1_correct = False
        read2_correct = False
        print("Error - check_pair - name of reads do not match: " + read1[0] + " : " + read2[0])

    if read1_correct and read2_correct:
        return True
    else:
        return False


def get_filename_stub(name, ext):
    stub = name

    pos = name.rfind(ext)
    if pos > 0:
        stub = name[:pos]

    return stub


def get_readname_stub(name):
    # https://en.wikipedia.org/wiki/FASTQ_format
    # @HWUSI-EAS100R:6:73:941:1973#0/1
    # @HWUSI-EAS100R:6:73:941:1973#0/2
    # @EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG
    # @EAS139:136:FC706VJ:2:2104:15343:197393 2:Y:18:ATCACG
    # @SRR001666.1 071112_SLXA-EAS1_s_7:5:1:817:345 length=36

    pos = 100000000
    ends = [name.find("/1"), name.find("/2"), name.find(" ")]
    # ends.append(name.find(" 1"))
    # ends.append(name.find(" 2"))

    for end in ends:
        if 0 < end < pos:
            pos = end

    stub = name

    if pos > 0:
        stub = name[:pos]

    return stub


def trim_name(name, delim):
    pos = len(name)
    if name.find(delim) > 0:
        pos = name.find(delim)

    return name[0:pos]


def check_fasta(line):
    if line.rfind(">") > 0:
        print("Error - FASTA sequence start symbol '>' found mid sequence " + line)


def check_seq_coding(seq, name):
    seq_test = False

    if len(seq) % 3 == 0:
        seq_test = True
    else:
        print("Sequence not divisible by 3: " + name)

    return seq_test


def fastq_to_fasta(filename):
    print("fastq_to_fasta")

    output_filename = filename.rfind(".") + ".fa"

    print("Input FASTQ file = " + filename)
    print("Output FASTA file =" + output_filename)

    seq_count = 0
    fastq_count = 0

    with open(output_filename, "w") as file_output:
        with open(filename) as file_handler:
            for line in file_handler:
                line = line.rstrip()
                fastq_count += 1

                if fastq_count == 1:
                    read = []

                read.append(line)

                if fastq_count == 4:
                    check_read(read)
                    read[0] = ">" + read[0][1:]
                    file_output.write(read[0] + "\n" + read[1] + "\n")
                    fastq_count = 0
                    seq_count += 1

        if fastq_count != 0:
            print("Error - fastq_to_fasta - incomplete FASTQ record in file (" + str(seq_count) + "/" + str(fastq_count) + "): " + filename)

    print("Sequences outputted = " + str(seq_count))


def fix_pairs(pair1_filename, pair2_filename):
    print("fix_pairs")

    pair1_out_filename = get_filename_stub(pair1_filename, ".f") + "_fixed.fq"
    pair2_out_filename = get_filename_stub(pair2_filename, ".f") + "_fixed.fq"
    singles1_out_filename = get_filename_stub(pair1_filename, ".f") + "_singles.fq"
    singles2_out_filename = get_filename_stub(pair2_filename, ".f") + "_singles.fq"

    print("Pair1 input file = " + pair1_filename)
    print("Pair2 input file = " + pair2_filename)
    print("Fixed Pair1 output file = " + pair1_out_filename)
    print("Fixed Pair2 output file = " + pair2_out_filename)
    print("Singletons1 output file = " + singles1_out_filename)
    print("Singletons2 output file = " + singles2_out_filename)

    seq_count = 0
    fastq_count = 0

    paired = 0
    single1 = 0
    single2 = 0

    reads1 = import_fastq(pair1_filename)

    with open(pair1_out_filename, "w") as pair1_out, open(pair2_out_filename, "w") as pair2_out, open(singles1_out_filename, "w") as singles1_out, open(singles2_out_filename, "w") as singles2_out:
        with open(pair2_filename) as file_handler:
            for line in file_handler:
                line = line.rstrip()
                fastq_count += 1

                if fastq_count == 1:
                    read2 = [get_readname_stub(line)]

                read2.append(line)

                if fastq_count == 4:
                    check_read(read2[1:])
                    found_pair = False

                    for read1 in reads1:
                        if read1[0] == read2[0]:
                            pair1_out.write(read1[1] + "\n" + read1[2] + "\n" + read1[3] + "\n" + read1[4] + "\n")
                            pair2_out.write(read2[1] + "\n" + read2[2] + "\n" + read2[3] + "\n" + read2[4] + "\n")
                            found_pair = True
                            reads1.remove(read1)
                            paired += 1
                            break

                    if not found_pair:
                        singles2_out.write(read2[1] + "\n" + read2[2] + "\n" + read2[3] + "\n" + read2[4] + "\n")
                        single2 += 1

                    fastq_count = 0
                    seq_count += 1

        for read1 in reads1:
            singles1_out.write(read1[1] + "\n" + read1[2] + "\n" + read1[3] + "\n" + read1[4] + "\n")
            single1 += 1

        print("Paired reads outputted = " + str(paired))
        print("Pair1-Singletons outputted = " + str(single1))
        print("Pair2-Singletons outputted = " + str(single2))


def fasta_add_length(filename):
    print("fasta_add_length")

    output_filename = get_filename_stub(filename, ".f") + "_lengths.fa"
    print("Input FASTA file = " + filename)
    print("Output FASTA file = " + output_filename)

    last_line = count_lines(filename)
    line_count = 0
    seq_count = 0

    with open(output_filename, "w") as file_output:
        with open(filename) as file_handler:

            seq = ''
            name = ''

            for line in file_handler:
                line = line.rstrip()
                line_count += 1

                if line.find(">") != 0:
                    seq += line.upper()

                if line.find(">") == 0 or line_count == last_line:
                    if seq_count > 0:
                        if len(seq) > 0:
                            file_output.write(name + "_len" + str(len(seq)) + "\n" + seq + "\n")
                        else:
                            print("Error - fasta_add_length - sequence length = 0: " + name)

                        name = ''
                        seq = ''

                    if line.find(">") == 0:
                        name = line
                        seq_count += 1

                check_fasta(line)

    print("Sequences outputted = " + str(seq_count))


def update_refseq_cai(seq, name, ref_codons_counts, ref_codons, ref_aa_counts):
    if check_seq_coding(seq, name):
        i = 0
        while i < len(seq):
            this_codon = seq[i] + seq[i + 1] + seq[i + 2]
            ns = this_codon.count("N")

            if ns == 0:
                ref_codons_counts[this_codon] += 1
                aa = ref_codons[this_codon]
                ref_aa_counts[aa] += 1

            i += 3


def fasta_calculate_cai(filename, allow_ns):
    print("fasta_calculate_cai")

    output_filename_table = get_filename_stub(filename, ".f") + "_cu_table.txt"
    output_filename_cai = get_filename_stub(filename, ".f") + "_cai.txt"
    print("Input FASTA file = " + filename)
    print("Output CU table file = " + output_filename_table)
    print("Output CAI data file = " + output_filename_cai)

    last_line = count_lines(filename)
    line_count = 0
    seq_count = 0
    proc_seq_count = 0

    ref_codons = create_codon_dict()
    ref_codons_counts = create_codon_counter()
    ref_codons_weights = create_codon_counter()
    ref_codons_freqs = create_codon_counter()
    ref_codons_thousand = create_codon_counter()
    ref_codons_rcsu = create_codon_counter()

    ref_aa_counts = create_aa_counter()
    ref_aa_syn = create_aa_syn_dict()
    ref_aa_max_syn = create_aa_counter()
    ref_aa_max_rcsu = create_aa_counter()

    with open(filename) as file_handler:
        seq = ''
        name = ''

        for line in file_handler:
            line = line.rstrip()
            line_count += 1

            if line.find(">") != 0:
                seq += line.upper()

            if line.find(">") == 0 or line_count == last_line:
                if seq_count > 0:
                    seq = correct_seq(seq)

                    if allow_ns == 1 or seq.count("N") == 0:
                        update_refseq_cai(seq, name, ref_codons_counts, ref_codons, ref_aa_counts)
                        proc_seq_count += 1

                    name = ''
                    seq = ''

                if line.find(">") == 0:
                    name = line
                    seq_count += 1

            check_fasta(line)

    print("Calculated ref codon table: processed sequences = " + str(proc_seq_count)+", total sequences = " + str(seq_count))

    calculate_ref_cai(ref_codons, ref_codons_counts, ref_codons_weights, ref_codons_thousand, ref_codons_freqs, ref_codons_rcsu, ref_aa_counts, ref_aa_syn, ref_aa_max_syn, ref_aa_max_rcsu)

    with open(output_filename_table, "w") as file_output:
        file_output.write("Codon\tCount\tWeight\tCodonBias\tNumberPer1000\tRCSU\n")

        for codon in ref_codons_counts:
            file_output.write(codon + "\t" + str(ref_codons_counts[codon]) + "\t" + str(ref_codons_weights[codon]) + "\t" + str(ref_codons_freqs[codon]) + "\t" + str(ref_codons_thousand[codon]) + "\t" + str(ref_codons_rcsu[codon]) + "\n")

        file_output.write("\n\n")
        file_output.write("Codon\tNumberPer1000\t(Count)\tCodon\tNumberPer1000\t(Count)\tCodon\tNumberPer1000\t(Count)\tCodon\tNumberPer1000\t(Count)\n")

        bases = ["T", "C", "A", "G"]

        for pos1 in bases:
            for pos3 in bases:
                counter = 0
                for pos2 in bases:
                    this_codon = pos1 + pos2 + pos3
                    uracil_codon = this_codon.replace("T", "U")

                    if counter > 0:
                        file_output.write("\t")

                    file_output.write(uracil_codon + "\t" + str(ref_codons_thousand[this_codon]) + "\t(" + str(ref_codons_counts[this_codon]) + ")")
                    counter += 1

                file_output.write("\n")

            file_output.write("\n")

    print("Outputted CU table")

    seq_count = 0
    proc_seq_count = 0
    line_count = 0

    with open(output_filename_cai, "w") as file_output:
        file_output.write("SeqName\tCAI\n")

        with open(filename) as file_handler:
            seq = ''
            name = ''

            for line in file_handler:
                line = line.rstrip()
                line_count += 1

                if line.find(">") != 0:
                    seq += line.upper()

                if line.find(">") == 0 or line_count == last_line:
                    if seq_count > 0:
                        seq = correct_seq(seq)
                        if allow_ns == 1 or seq.count("N") == 0:
                            this_cai = calculate_seq_cai(seq, name, ref_codons, ref_codons_weights)
                            file_output.write(name + "\t" + str(this_cai) + "\n")
                            proc_seq_count += 1

                        name = ''
                        seq = ''

                    if line.find(">") == 0:
                        name = line
                        seq_count += 1

                check_fasta(line)

        print("Calculated CAI for each sequence: processed sequences = " + str(proc_seq_count)+", total sequences = " + str(seq_count))


def clean_seq_cai(seq, seq_stats):
    # seq_stats: start, stop, mid_stop, ncodon, incomplete, short

    ref_codons = create_codon_dict()

    seq = correct_seq(seq)

    new_seq = ""

    if len(seq) >= 3:
        this_codon = seq[0:3]

        if this_codon == "ATG":
            seq_stats[0] += 1

        i = 0
        while i <= len(seq) - 3:
            this_codon = seq[i] + seq[i + 1] + seq[i + 2]
            new_seq += this_codon
            i += 3

        if i < len(seq):
            seq_stats[4] += 1
    else:
        new_seq = ""

    seq = new_seq
    new_seq = ""

    last_codon = seq[-3:]
    if ref_codons[last_codon] == "X":
        seq = seq[:-3]
        seq_stats[1] += 1

    if len(seq) >= 3:
        i = 0
        while i < len(seq):
            this_codon = seq[i] + seq[i + 1] + seq[i + 2]

            if this_codon.count("N") == 0:
                if ref_codons[this_codon] == "X":
                    new_seq += "NNN"
                    seq_stats[2] += 1
                else:
                    new_seq += this_codon
            else:
                new_seq += "NNN"
                seq_stats[3] += 1

            i += 3
    else:
        new_seq = ""

    seq = new_seq

    if len(seq) < 3:
        seq_stats[5] += 1

    return seq


def calculate_seq_cai(seq, name, ref_codons, ref_codons_weights):
    weight = 0
    this_codon_count = 0

    # check divisible by 3
    if check_seq_coding(seq, name):
        i = 0
        while i < len(seq):
            this_codon = seq[i] + seq[i + 1] + seq[i + 2]

            ns = this_codon.count("N")

            if ns == 0:
                this_aa = ref_codons[this_codon]
                if this_aa != "M" and this_aa != "W" and this_aa != "X":
                    weight += log(ref_codons_weights[this_codon])
                    this_codon_count += 1

            i += 3

    if this_codon_count == 0:
        print("Warning - calculate_seq_cai - sequence has 0 usable codons: " + name + " " + seq)
    else:
        weight = exp(weight / this_codon_count)

    return weight


def calculate_ref_cai(ref_codons, ref_codons_counts, ref_codons_weights, ref_codons_thousand, ref_codons_freqs, ref_codons_rcsu, ref_aa_counts, ref_aa_syn, ref_aa_max_syn, ref_aa_max_rcsu):
    total_codons = 0

    for codon in ref_codons_counts:
        total_codons += ref_codons_counts[codon]

    for codon in ref_codons_counts:
        ref_codons_thousand[codon] = ref_codons_counts[codon] / total_codons * 1000

    for aa in ref_aa_max_syn:
        this_max = 0

        for codon in ref_codons:
            if ref_codons[codon] == aa:
                if ref_codons_counts[codon] > this_max:
                    this_max = ref_codons_counts[codon]

        ref_aa_max_syn[aa] = this_max

    for aa in ref_aa_max_syn:
        for codon in ref_codons:
            if ref_codons[codon] == aa:
                if ref_codons_counts[codon] > 0:
                    ref_codons_weights[codon] = ref_codons_counts[codon] / ref_aa_max_syn[aa]
                else:
                    ref_codons_weights[codon] = 0

    for codon in ref_codons_freqs:
        if ref_aa_counts[ref_codons[codon]] > 0:
            ref_codons_freqs[codon] = ref_codons_counts[codon] / ref_aa_counts[ref_codons[codon]]

        elif ref_codons_counts[codon] > 0:
            print("Error - calculate_ref_cai - codon count > 0 but AA count = 0: "+codon)

    for codon in ref_codons_counts:
        aa = ref_codons[codon]
        if ref_aa_counts[aa] > 0:
            ref_codons_rcsu[codon] = ref_codons_counts[codon] / (ref_aa_counts[aa] * (1 / ref_aa_syn[aa]))

    for aa in ref_aa_max_rcsu:
        this_max = 0

        for codon in ref_codons:
            if ref_codons[codon] == aa:
                if ref_codons_rcsu[codon] > this_max:
                    this_max = ref_codons_rcsu[codon]

        ref_aa_max_rcsu = this_max


def fasta_filter_length(filename, min_len):
    print("fasta_filter_length")

    output_filename = get_filename_stub(filename, ".f") + "_" + str(min_len) + ".fa"
    print("Input FASTA file = " + filename)
    print("Minimum sequence length = " + str(min_len))
    print("Output FASTA file = " + output_filename)

    last_line = count_lines(filename)
    line_count = 0
    seq_count = 0
    out_count = 0
    not_count = 0

    with open(output_filename, "w") as file_output:
        with open(filename) as file_handler:

            seq = ''
            name = ''

            for line in file_handler:
                line = line.rstrip()
                line_count += 1

                if line.find(">") != 0:
                    seq += line.upper()

                if line.find(">") == 0 or line_count == last_line:
                    if seq_count > 0:
                        if len(seq) >= min_len:
                            file_output.write(name + "\n" + seq + "\n")
                            out_count += 1
                        else:
                            not_count += 1

                        name = ''
                        seq = ''

                    if line.find(">") == 0:
                        name = line
                        seq_count += 1

                check_fasta(line)

    print("Sequences outputted >= " + str(min_len) + " = " + str(out_count))
    print("Sequences < " + str(min_len) + " = " + str(not_count))


def fasta_rna_to_dna(filename):
    print("fasta_rna_to_dna")

    output_filename = get_filename_stub(filename, ".f") + "_dna.fa"
    print("Input RNA FASTA file = " + filename)
    print("Output DNA FASTA file = " + output_filename)

    seq_count = 0
    line_count = 0

    with open(output_filename, "w") as file_output:
        with open(filename) as file_handler:
            for line in file_handler:
                line = line.rstrip()
                line_count += 1

                if line.find(">") == 0:
                    file_output.write(line + "\n")
                    seq_count += 1
                else:
                    file_output.write(rna_to_dna(line) + "\n")

                check_fasta(line)

    print("Sequences outputted = " + str(seq_count))


def fasta_dna_to_rna(filename):
    print("fasta_dna_to_rna")

    output_filename = get_filename_stub(filename, ".f") + "_rna.fa"
    print("Input DNA FASTA file = " + filename)
    print("Output RNA FASTA file = " + output_filename)

    seq_count = 0
    line_count = 0

    with open(output_filename, "w") as file_output:
        with open(filename) as file_handler:
            for line in file_handler:
                line = line.rstrip()
                line_count += 1

                if line.find(">") == 0:
                    file_output.write(line + "\n")
                    seq_count += 1
                else:
                    file_output.write(dna_to_rna(line) + "\n")

                check_fasta(line)

    print("Sequences outputted = " + str(seq_count))


def fasta_trim_names(filename, delim=" "):
    print("fasta_trim_names")

    output_filename = get_filename_stub(filename, ".f") + "_trim.fa"
    print("Input FASTA file = " + filename)
    print("Delimiter = '" + delim + "'")
    print("Output FASTA file = " + output_filename)

    seq_count = 0
    line_count = 0

    with open(output_filename, "w") as file_output:
        with open(filename) as file_handler:
            for line in file_handler:
                line = line.rstrip()
                line_count += 1

                if line.find(">") == 0:
                    file_output.write(trim_name(line, delim) + "\n")
                    seq_count += 1
                else:
                    file_output.write(line + "\n")

                check_fasta(line)

    print("Sequences outputted = " + str(seq_count))


def subsample_pairs(pair1_filename, pair2_filename, sub_reads):
    print("subsample_pairs")

    pair1_out_filename = get_filename_stub(pair1_filename, ".f") + "_N" + str(sub_reads) + ".fq"
    pair2_out_filename = get_filename_stub(pair2_filename, ".f") + "_N" + str(sub_reads) + ".fq"

    print("Pair1 input file = " + pair1_filename)
    print("Pair2 input file = " + pair2_filename)
    print("Number of reads to subsample = " + sub_reads)
    print("Pair1 output file = " + pair1_out_filename)
    print("Pair2 output file = " + pair2_out_filename)

    num_reads_1 = count_fastq(pair1_filename)
    num_reads_2 = count_fastq(pair1_filename)

    print("Number of read in Pair1 = " + str(num_reads_1))
    print("Number of read in Pair2 = " + str(num_reads_2))

    if num_reads_1 != num_reads_2:
        print("Error - number or reads in pair1 and pair2 not equal")
        sys.exit(1)

    if sub_reads >= num_reads_1:
        print("Error - subsample_pairs - number or reads to be subsampled (" + str(sub_reads) + ") >= the number of reads (" + str(num_reads_1))
        sys.exit(1)

    all_reads = range(0, num_reads_1)
    sel_reads = random.sample(all_reads, sub_reads)

    seq_count = 0
    fastq_count = 0
    pair_out_count = 0
    pair_not_count = 0
    pair_error_count = 0

    with open(pair1_out_filename, "w") as pair1_out, open(pair2_out_filename, "w") as pair2_out:
        with open(pair1_filename) as pair1_handler, open(pair2_filename) as pair2_handler:
            for line1, line2 in zip(pair1_handler, pair2_handler):
                line1 = line1.rstrip()
                line2 = line2.rstrip()

                fastq_count += 1

                if fastq_count == 1:
                    read1 = [get_readname_stub(line1)]
                    read2 = [get_readname_stub(line2)]

                read1.append(line1)
                read2.append(line2)

                if fastq_count == 4:
                    if seq_count in sel_reads:
                        if check_pair(read1, read2):
                            pair1_out.write(read1[1] + "\n" + read1[2] + "\n" + read1[3] + "\n" + read1[4] + "\n")
                            pair2_out.write(read2[1] + "\n" + read2[2] + "\n" + read2[3] + "\n" + read2[4] + "\n")
                            pair_out_count += 1
                        else:
                            print("Error - subsample_pairs - could not output the " + str(seq_count) + " read pair as read names not equal")
                            pair_error_count += 1
                    else:
                        pair_not_count += 1

                    fastq_count = 0
                    seq_count += 1

        print("Pairs outputted = " + str(pair_out_count))
        print("Erroneous pairs not outputted = " + str(pair_error_count))


def subsample_reads(read_filename, sub_reads):
    print("subsample_reads")

    read_out_filename = get_filename_stub(read_filename, ".f") + "_N" + str(sub_reads) + ".fq"
    print("Input FASTQ file = " + read_filename)
    print("Number of reads to subsample = " + str(sub_reads))
    print("Output FASTQ file = " + read_out_filename)

    num_reads = count_fastq(read_filename)

    if sub_reads >= num_reads:
        print("Error - subsample_reads - number or reads to be subsampled (" + str(sub_reads) + ") >= the number of reads (" + str(num_reads))
        sys.exit(1)

    print("Number of reads = " + str(num_reads))

    all_reads = range(0, num_reads)
    sel_reads = random.sample(all_reads, sub_reads)

    seq_count = 0
    fastq_count = 0
    read_out_count = 0
    read_not_count = 0
    read_error_count = 0

    with open(read_out_filename, "w") as read_out:
        with open(read_filename) as read_handler:
            read = []

            for line in read_handler:
                line = line.rstrip()

                fastq_count += 1

                if fastq_count == 1:
                    read = []

                if 1 <= fastq_count <= 4:
                    read.append(line)

                if fastq_count == 4:
                    if seq_count in sel_reads:
                        if check_read(read):
                            read_out.write(read[0] + "\n" + read[1] + "\n" + read[2] + "\n" + read[3] + "\n")
                            read_out_count += 1
                        else:
                            print("Error - subsample_reads - could not output the " + str(seq_count) + " read as not in correct FASTQ format")
                            read_error_count += 1
                    else:
                        read_not_count += 1

                    fastq_count = 0
                    seq_count += 1

        print("Outputted read pairs = " + str(read_out_count))
        print("Erroneous pairs not outputted = " + str(read_error_count))


def fasta_select_seqs(id_filename, id_column, seq_filename):

    print("fasta_select_seqs")
    print("Sequence ID file = " + id_filename)
    print("Sequence ID column = " + str(id_column))
    print("Input sequence file = " + seq_filename)

    output_filename = get_filename_stub(seq_filename, ".f") + "_sel.fa"
    print("Output sequence file = " + output_filename)

    seq_ids = []
    seq_count = 0
    with open(id_filename) as id_handler:
        reader = csv.reader(id_handler, delimiter='\t')
        for row in reader:
            if seq_count > 0:
                seq_ids.append(row[id_column])

            seq_count += 1

    # header
    seq_count -= 1

    print("Sequences in ID file = " + str(seq_count))

    last_line = count_lines(seq_filename)

    line_count = 0
    seq_count = 0
    seq_found = 0

    # start, stop, mid_stop, ncodon, incomplete, short
    seq_stats = [0, 0, 0, 0, 0, 0]

    with open(output_filename, "w") as file_output:
        with open(seq_filename) as seq_handler:
            seq = ''
            name = ''

            for line in seq_handler:
                line = line.rstrip()
                line_count += 1

                if line.find(">") != 0:
                    seq += line.upper()

                if line.find(">") == 0 or line_count == last_line:
                    if seq_count > 0:
                        if name in seq_ids:
                            seq = clean_seq_cai(seq, seq_stats)
                            file_output.write(name + "\n" + seq + "\n")
                            seq_found += 1

                        name = ''
                        seq = ''

                    if line.find(">") == 0:
                        name = line
                        seq_count += 1

                check_fasta(line)

    print("Found seqs = " + str(seq_found))
    print("Seq with ATG at start = " + str(seq_stats[0]))
    print("Seq with STOP at end = " + str(seq_stats[1]))
    print("Number of mid-STOP codons = " + str(seq_stats[2]))
    print("Number of N codons = " + str(seq_stats[3]))
    print("Number of incomplete seqs = " + str(seq_stats[4]))
    print("Number of seqs <3 = " + str(seq_stats[5]))


def fasta_remove_duplicates(filename):
    print("fasta_remove_duplicates")

    output_filename = get_filename_stub(filename, ".f") + "_nodups.fa"
    print("Input FASTA file = " + filename)
    print("Output FASTA file = " + output_filename)

    last_line = count_lines(filename)
    line_count = 0
    seq_count = 0
    out_count = 0
    dup_count = 0

    all_seqs = {}

    with open(filename) as file_handler:

        seq = ''
        name = ''

        for line in file_handler:
            line = line.rstrip()
            line_count += 1

            if line.find(">") != 0:
                seq += line.upper()

            if line.find(">") == 0 or line_count == last_line:
                if seq_count > 0:
                    if name in all_seqs:
                        dup_count += 1
                        if all_seqs[name] != seq:
                            print("Error - duplicate seq names but different seqs: " + name + " - " + seq)
                    else:
                        all_seqs[name] = seq

                    name = ''
                    seq = ''

                if line.find(">") == 0:
                    name = line
                    seq_count += 1

            check_fasta(line)

    with open(output_filename, "w") as file_output:
        for seq in all_seqs:
            file_output.write(seq+"\n"+all_seqs[seq]+"\n")
            out_count += 1

    print("Sequences = " + str(seq_count))
    print("Duplicates = " + str(dup_count))
    print("Outputted = " + str(out_count))


def fasta_length_distribution(filename):
    print("fasta_length_distribution")

    print("Input FASTA file = " + filename)
    last_line = count_lines(filename)
    line_count = 0
    seq_count = 0

    seq_lengths = {300000: 0, 250000: 0, 200000: 0, 150000: 0, 100000: 0, 50000: 0, 25000: 0, 10000: 0, 5000: 0, 1000: 0, 0: 0}

    with open(filename) as file_handler:
        seq = ''

        for line in file_handler:
            line = line.rstrip()
            line_count += 1

            if line.find(">") != 0:
                seq += line.upper()

            if line.find(">") == 0 or line_count == last_line:
                if seq_count > 0:
                    for slen in seq_lengths:
                        if len(seq) >= slen:
                            seq_lengths[slen] += 1

                    seq = ''

                if line.find(">") == 0:
                    seq_count += 1

            check_fasta(line)

    print("Lenght distribution:")

    for slen in seq_lengths:
        print(str(slen) + "\t" + str(seq_lengths[slen]))


def fasta_align_distance(filename):
    print("fasta_align_distance")

    output_filename = get_filename_stub(filename, ".f") + "_dist.txt"
    print("Input FASTA file = " + filename)
    print("Output FASTA file = " + output_filename)

    last_line = count_lines(filename)
    line_count = 0
    seq_count = 0

    all_seqs = []

    with open(filename) as file_handler:
        seq = ''
        name = ''

        for line in file_handler:
            line = line.rstrip()
            line_count += 1

            if line.find(">") != 0:
                seq += line.upper()

            if line.find(">") == 0 or line_count == last_line:
                if seq_count > 0:
                    all_seqs.append([name, seq])

                    name = ''
                    seq = ''

                if line.find(">") == 0:
                    name = line
                    seq_count += 1

            check_fasta(line)

    with open(output_filename, "w") as file_out:
        for seq in all_seqs:
            file_out.write("\t" + seq[0])
        file_out.write("\n")

        for seq in all_seqs:
            if len(seq[1]) != len(all_seqs[0][1]):
                print("Error - seq in alignment not the same length: " + seq[0])
            else:
                file_out.write(seq[0])

                for seq2 in all_seqs:
                    seq_dist = 0
                    seq_dist_nogaps = 0

                    for base in range(len(seq2[1])):
                        if seq2[1][base] != seq[1][base]:
                            seq_dist += 1

                            if base != "-" and seq[1][base] != "-":
                                seq_dist_nogaps += 1
                    # probably need a length without GAPs as well

                    seq_dist = seq_dist / len(all_seqs[0][1])
                    seq_dist_nogaps = seq_dist_nogaps / len(all_seqs[0][1])

                    dist = "%.8f" % seq_dist
                    file_out.write("\t" + dist)

            file_out.write("\n")


print("ortools.py started...")

example_messages = ["ortools.py fix_pairs pair1_file pair2_file",
                    "ortools.py subsample_pairs pair1_file pair2_file number_of_reads",
                    "ortools.py subsample_reads read_file number_of_reads",
                    "ortools.py fasta_filter_length fasta_file min_len",
                    "ortools.py fasta_add_length fasta_file",
                    "ortools.py fasta_rna_to_dna fasta_file",
                    "ortools.py fasta_dna_to_rna fasta_file",
                    "ortools.py fasta_trim_names fasta_file",
                    "ortools.py fasta_calculate_cai fasta_file",
                    "ortools.py fasta_select_seqs id_file id_column seq_file",
                    "ortools.py fasta_remove_duplicates fasta_file",
                    "ortools.py fasta_length_distribution fasta_file",
                    "ortools.py fasta_align_distance fasta_file",
                    "ortools.py fastq_to_fasta fastq_file"]

error_message = "Error - incorrect number of arguments - example usage:"

arguments = len(sys.argv)

if arguments == 1:
    for example in example_messages:
        print(example)

elif 1 < arguments:

    if sys.argv[1] == "fix_pairs":
        if arguments != 4:
            print(error_message + "\n" + example_messages[0])
        else:
            fix_pairs(sys.argv[2], sys.argv[3])

    elif sys.argv[1] == "subsample_pairs":
        if arguments != 5:
            print(error_message + "\n" + example_messages[1])
        else:
            try:
                arg_number_of_reads = int(sys.argv[4])
            except ValueError:
                print("Error - unrecognised integer number for number_of_reads to sample: "+sys.argv[4])
                sys.exit(1)

            subsample_pairs(sys.argv[2], sys.argv[3], arg_number_of_reads)

    elif sys.argv[1] == "subsample_reads":
        if arguments != 4:
            print(error_message + "\n" + example_messages[2])
        else:
            try:
                arg_number_of_reads = int(sys.argv[3])
            except ValueError:
                print("Error - unrecognised integer number for number_of_reads to sample: "+sys.argv[3])
                sys.exit(1)

            subsample_reads(sys.argv[2], arg_number_of_reads)

    elif sys.argv[1] == "fasta_filter_length":
        if arguments != 4:
            print(error_message + "\n" + example_messages[3])
        else:
            try:
                arg_len = int(sys.argv[3])
            except ValueError:
                print("Error - unrecognised integer number for length: "+sys.argv[3])
                sys.exit(1)

            fasta_filter_length(sys.argv[2], arg_len)

    elif sys.argv[1] == "fasta_add_length":
        if arguments != 3:
            print(error_message + "\n" + example_messages[4])
        else:
            fasta_add_length(sys.argv[2])

    elif sys.argv[1] == "fasta_rna_to_dna":
        if arguments != 3:
            print(error_message + "\n" + example_messages[5])
        else:
            fasta_rna_to_dna(sys.argv[2])

    elif sys.argv[1] == "fasta_dna_to_rna":
        if arguments != 3:
            print(error_message + "\n" + example_messages[6])
        else:
            fasta_dna_to_rna(sys.argv[2])

    elif sys.argv[1] == "fasta_trim_names":
        if arguments != 3:
            print(error_message + "\n" + example_messages[7])
        else:
            fasta_trim_names(sys.argv[2])

    elif sys.argv[1] == "fasta_calculate_cai":
        if arguments != 4:
            print(error_message + "\n" + example_messages[8])
        else:
            try:
                arg_ns = int(sys.argv[3])
            except ValueError:
                print("Error - unrecognised integer number for length: "+sys.argv[3])
                sys.exit(1)

            fasta_calculate_cai(sys.argv[2], arg_ns)

    elif sys.argv[1] == "fasta_select_seqs":
        if arguments != 5:
            print(error_message + "\n" + example_messages[9])
        else:
            try:
                arg_col = int(sys.argv[3])
            except ValueError:
                print("Error - unrecognised integer number for length: "+sys.argv[3])
                sys.exit(1)

            fasta_select_seqs(sys.argv[2], arg_col, sys.argv[4])

    elif sys.argv[1] == "fasta_remove_duplicates":
        if arguments != 3:
            print(error_message + "\n" + example_messages[10])
        else:
            fasta_remove_duplicates(sys.argv[2])

    elif sys.argv[1] == "fasta_length_distribution":
        if arguments != 3:
            print(error_message + "\n" + example_messages[11])
        else:
            fasta_length_distribution(sys.argv[2])

    elif sys.argv[1] == "fasta_align_distance":
        if arguments != 3:
            print(error_message + "\n" + example_messages[12])
        else:
            fasta_align_distance(sys.argv[2])

    elif sys.argv[1] == "fastq_to_fasta":
        if arguments != 3:
            print(error_message + "\n" + example_messages[13])
        else:
            fastq_to_fasta(sys.argv[2])

# convert fastatofastq [default q40]
# filter fastq by length

# remove sequence duplicates (based on full name, or up to a delimiter)
# remove duplicates from text file - specify column to use a duplicate

# AlignTools
# AlignMuts
# Consenus of Alignment
# Remove any site that has a GAP
# Hamming distance
# Coding sequence

# VCFCOnverter
# diversi, varscan, vphaser

# https://biopython.org/DIST/docs/api/Bio.SeqIO-module.html
# https://biopython.org/wiki/SeqIO
# https://thispointer.com/5-different-ways-to-read-a-file-line-by-line-in-python/

# addLengthToSequenceName("/Users/richardorton/Downloads/contigs.fasta")
# fix_fastq_pairs("/Users/richardorton/Downloads/404_1406776_S90_R1_001.fastq",
#                "/Users/richardorton/Downloads/404_1406776_S90_R2_001.fastq")
# removeFastaSequencesByLength("/Users/richardorton/Downloads/contigs.fasta",1000)

# sys.stdout.write("Hey! Welcome to the STechies.\n")
# print("Good Morning!", end = '')

# String iteration
# https://thispointer.com/python-how-to-iterate-over-the-characters-in-string/

print("...finished ortools.py")
