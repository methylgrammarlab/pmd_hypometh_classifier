import os
import re

# output/input file format
BULK_FILE_FORMAT = "all_cpg_ratios_%s_hg19.dummy.pkl.zip"
SCWGBS_FILE_FORMAT = "all_cpg_ratios_%s_%s_hg19.dummy.pkl.zip"
SCWGBS_FILE_FORMAT = "all_cpg_ratios_%s_%s.dummy.pkl.zip"  # TODO :use the previous one
FULL_CPG_CONTEXT_FILE_FORMAT = "full_cpg_seq_chr*_hg19.pickle.zlib"

# regex
DATA_FILE_SCWGBS_RE = re.compile(".+(CRC\d+)_(chr\d+).dummy.pkl.zip")
DATA_FILE_BULK_RE = re.compile(".+_(chr\d+)_hg19.dummy.pkl.zip")

PATIENT_CHR_NAME_RE = re.compile(".*(CRC\d+)_(chr\d+).*")

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
LABEL_HYPO_RESISTANCE = 0  # Group 0 - partial loss
LABEL_HYPO_PRONE = 1  # Group 1 - comp loss

# Parsing options
SCWGBS = "scwgbs"
SCWGBS_CRC01 = "scwgbs_CRC01"
BULK = "bulk"
PARSE_FORMAT_OPTIONS = [SCWGBS, BULK,SCWGBS_CRC01]

# IMPORTANT FILES
MAIN_FOLDER = r"/vol/sci/bio/data/benjamin.berman/bermanb/projects/scTrio-seq-reanalysis/liordror/"
CONTEXT_MAP_FULL = os.path.join(MAIN_FOLDER, r"genomic_data/full")
CONTEXT_MAP_FILTERED_NO_BL_CPGI = os.path.join(MAIN_FOLDER, r"/genomic_data/filtered_by_bl_and_cpgi")
SUBLINEAGE_FILE = os.path.join(MAIN_FOLDER, r"/genomic_data/orig/sublineage/patient_sublineage_dict.pickle")
PMD_FILE = os.path.join(MAIN_FOLDER, r"genomic_data/orig/pmds/pmd_dict.pickle")
GENOME_FILE = os.path.join(MAIN_FOLDER, "/genomic_data/orig/seq_context/genome.fa")
NC_FILES = os.path.join(MAIN_FOLDER, "stats/average_nc_methylation/all")



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
CONVERT_SUBLINEAGE_LIOR = R"C:\Users\liorf\OneDrive\Documents\University\year 3\Project\proj_scwgbs\samples_handling\my_files\top_bottom\convert_sublineage.pickle.zlib"
CONVERT_SUBLINEAGE_LIOR_AQUARIUM = R"/cs/usr/liorf/PycharmProjects/proj_scwgbs/stats/convert_sublineage.pickle.zlib"
NC_LIOR = r"C:\Users\liorf\Documents\average_nc_methylation"
GENOME_FILE_LOCAL_LIOR = r"C:\Users\liorf\Documents\genome.fa"
