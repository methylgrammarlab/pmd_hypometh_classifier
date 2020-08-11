"""
Several functions related to sequence handling
"""
import tqdm

from commons import consts

TRANSLATION_TABLE = {84: 65, 65: 84, 67: 71, 71: 67}
SEQ_SIZE = 150

try:
    import pyfaidx

    genome = pyfaidx.Fasta(filename=consts.GENOME_FILE_LOCAL_DROR, sequence_always_upper=True, as_raw=True)
except ImportError:
    print("Error: no pyfaidx install, somme code might fail")


def get_reverse_seq(sequences):
    """
    Return a list of reverse strand sequences
    :param sequences: some kind of iter with the strands to reverse
    :return: A list with the reverse strand
    """
    return [s.translate(TRANSLATION_TABLE)[::-1] for s in sequences]


def get_cpg_sequence(chr_num, cpg_index, seq_size=SEQ_SIZE):
    """
    Get the sequence of a specific length for a CpG
    :param chr_num: The name of the chromosome, needs to be "16" and not int
    :type chr_num: str
    :param cpg_index: The cpg index
    :param seq_size: The sequence size
    :return: The sequence with the required length with the CG in the middle
    """
    chr_info = genome[chr_num]
    seq_size = int((seq_size - 2) / 2)
    return chr_info[cpg_index - seq_size - 1:cpg_index + seq_size + 1]


def get_sequences_for_cpgs(indexes_list, chromosome, seq_size=SEQ_SIZE):
    """
    Get all the sequences for the required cpg
    :param indexes_list: An list of indexes for cpg
    :param chromosome: The chromosome for those cpg
    :param seq_size: The sequence size
    :return: A list with all the sequences
    """
    if isinstance(chromosome, int):
        chromosome = str(chromosome)
    elif isinstance(chromosome, str):
        if chromosome.startswith("chr"):
            chromosome = chromosome[3:]

    sequences = []
    for index in tqdm.tqdm(indexes_list):
        sequences.append(get_cpg_sequence(chr_num=chromosome, cpg_index=index, seq_size=seq_size))

    return sequences
