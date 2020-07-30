import os
import re

# output/input file format
BULK_FILE_FORMAT = "all_cpg_ratios_%s_hg19.dummy.pkl.zip"
SCWGBS_FILE_FORMAT = "all_cpg_ratios_%s_%s.dummy.pkl.zip"
FULL_CPG_CONTEXT_FILE_FORMAT = "full_cpg_seq_chr*.pickle.zlib"

# regex
DATA_FILE_SCWGBS_RE = re.compile(".+(CRC\d+)_(chr\d+)_hg19.dummy.pkl.zip")
DATA_FILE_BULK_RE = re.compile(".+_(chr\d+)_hg19.dummy.pkl.zip")

PATIENT_CHR_NUM_RE = re.compile(".*(CRC\d+)_chr(\d+).*")
PATIENT_CHR_NAME_RE = re.compile(".*(CRC\d+)_(chr\d+).*")

CHR_NUM_FROM_FULL_RE = re.compile("chr(\d+).*")
CHR_FULL_NAME_RE = re.compile("(chr\d+)")

# Regex for weak and strong context cpg
WEAK_FORMAT_RE = re.compile("\w\w[AT]CG[AT]\w\w")
STRONG_FORMAT_RE = re.compile("\w\w[CG]CG[CG]\w\w")

# Mapping between nucleotide and number
NUCLEOTIDE_TO_NUMBER = {"A": "1",
                        "C": "2",
                        "G": "3",
                        "T": "4",
                        "N": "5"}  # N means that this is unknown

# NN labels
LABEL_PARTIAL_LOST = 0  # Group 0 - partial loss
LABEL_COMPLETELY_LOST = 1  # Group 1 - comp loss


CENETROMERE_DICT = {
    "1": [121311574, 129523877],
    "2": [91278726, 97911435],
    "3": [88724336, 93610603],
    "4": [48396285, 52374061],
    "5": [45111338, 50515299],
    "6": [57556887, 63334797],
    "7": [57041911, 63035444],
    "8": [43719124, 48091035],
    "9": [47315670, 51717126],
    "10": [38196156, 42596634],
    "11": [52073942, 56281937],
    "12": [33028390, 38938733],
    "13": [16004126, 20341692],
    "14": [16172139, 19796928],
    "15": [15579446, 20905751],
    "16": [34499088, 39310184],
    "17": [22776839, 25729391],
    "18": [15615450, 19164415],
    "19": [24342712, 28719791],
    "20": [25701316, 29630179],
    "21": [10688588, 14563981],
    "22": [12326422, 17790024],
}

# IMPORTANT FILES
MAIN_FOLDER = r"/vol/sci/bio/data/benjamin.berman/bermanb/projects/scTrio-seq-reanalysis/liordror/"
CONTEXT_MAP_FULL = os.path.join(MAIN_FOLDER, r"genomic_data/full")
CONTEXT_MAP_FILTERED_NO_BL_CPGI = os.path.join(MAIN_FOLDER, r"/genomic_data/filtered_by_bl_and_cpgi")
SUBLINEAGE_FILE = os.path.join(MAIN_FOLDER, r"/genomic_data/orig/sublineage/patient_sublineage_dict.pickle")
PMD_FILE = os.path.join(MAIN_FOLDER, r"genomic_data/orig/pmds/pmd_dict.pickle")
GENOME_FILE = os.path.join(MAIN_FOLDER, "/genomic_data/orig/seq_context/genome.fa")



# DROR local
CONTEXT_MAP_FILTERED_LOCAL_DROR = r"H:\Study\university\Computational-Biology\Year 3\Projects\proj_scwgbs\resource\genomic_data\filtered_by_bl_and_cpgi"
PMD_FILE_LOCAL_DROR = r"H:\Study\university\Computational-Biology\Year 3\Projects\proj_scwgbs\resource\pmd_dict.pickle"
SUBLINEAGE_FILE_LOCAL_DROR = r"H:\Study\university\Computational-Biology\Year " \
                             r"3\Projects\proj_scwgbs\resource\sublineage\patient_sublineage_dict.pickle"
GENOME_FILE_LOCAL_DROR = r"H:\Study\university\Computational-Biology\Year " \
                         r"3\Projects\proj_scwgbs\resource\full_genome\genome.fa"
CONVERT_SUBLINEAGE_DROR = R"H:\Study\university\Computational-Biology\Year 3\Projects\proj_scwgbs\resource\convert_sublineage.pickle.zlib"
NC_DROR = r"H:\Study\university\Computational-Biology\Year 3\Projects\proj_scwgbs\resource\average_nc_methylation"
CONTEXT_MAP_FILTERED_NO_BL_CPGI_DROR = r"H:\Study\university\Computational-Biology\Year " \
                                       r"3\Projects\proj_scwgbs\resource\genomic_data\all"

# Lior local
PMD_FILE_LOCAL_LIOR = r"C:\Users\liorf\OneDrive\Documents\University\year 3\Project\proj_scwgbs\resources\pmd_dict.pickle"
CONTEXT_MAP_FILTERED_NO_BL_CPGI_LIOR = r"C:\Users\liorf\Documents\genomic_data\filtered_by_bl_and_cpgi"
CONVERT_SUBLINEAGE_LIOR = R"C:\Users\liorf\OneDrive\Documents\University\year 3\Project\proj_scwgbs\stats\top_bottom\convert_sublineage.pickle.zlib"
CONVERT_SUBLINEAGE_LIOR_AQUARIUM = R"/cs/usr/liorf/PycharmProjects/proj_scwgbs/stats/convert_sublineage.pickle.zlib"

