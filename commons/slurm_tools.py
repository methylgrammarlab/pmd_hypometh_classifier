import datetime
import glob
import os
import re
import sys

CPG_FORMAT_FILE_RE = re.compile(".+(CRC\d+)_(chr\d+).dummy.pkl.zip")


def init_slurm(func):
    """
    Some simple prints when running code in slurm
    :param func:
    :return:
    """
    print("Run the script: %s" % " ".join(sys.argv))
    start_time = datetime.datetime.now()
    print("Start time: %s" % start_time)
    try:
        func()
    except Exception as ex:
        print("Run ended due to exception: %s" % ex)
    else:
        print("Run ended with no exceptions")
    end_time = datetime.datetime.now()
    delta = end_time - start_time
    print("End time:%s\nTotal time is:%s seconds" % (end_time, delta.total_seconds()))


def create_sbatch():
    code = """#!/bin/bash
    #SBATCH --output=/cs/usr/drordod/Desktop/project/proj_scwgbs/output/coverage_between_pairs_%s.out
    #SBATCH --error=/cs/usr/drordod/Desktop/project/proj_scwgbs/output/coverage_between_pairs_%s.err
    #SBATCH --time=72:00:00
    #SBATCH --mem=4G
    
    RUNPATH=~/Desktop/project/proj_scwgbs
    cd $RUNPATH/format_files
    python3 coverage_between_pairs.py --cpg_format_files %s --output_folder /cs/usr/drordod/Desktop/scTrio-seq-reanalysis/liordror/stats/coverage_histograms/
    """

    input_folder = "/cs/usr/drordod/Desktop/scTrio-seq-reanalysis/liordror/cpg_format/filtered"
    all_files = glob.glob(os.path.join(input_folder, "*", "all_cpg_ratios_*.dummy.pkl.zip"))
    output_path = r"/cs/usr/drordod/Desktop/project/proj_scwgbs/sbatch/coverage_between_pairs"
    names = []

    for cpg_format_file in all_files:
        patient, chromosome = CPG_FORMAT_FILE_RE.findall(cpg_format_file)[0]
        label = "%s_%s" % (patient, chromosome)
        names.append("%s.sbatch" % label)
        output_f = os.path.join(output_path, "%s.sbatch" % (label))
        with open(output_f, "w") as output:
            data = code % (label, label, cpg_format_file)
            output.write(data)

    bin_path = os.path.join(output_path, "run_sbatch")
    with open(bin_path, "w") as bin_file:
        bin_file.write("#!/bin/bash\n")
        for name in names:
            bin_file.write("sbatch %s\n" % name)
