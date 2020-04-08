import pandas as pd

crc01_path = r"H:\Study\university\Computational-Biology\Year " \
             r"3\Projects\proj_scwgbs\resource\cpg_format\filtered_by_bl_and_cpgi\CRC01\all_cpg_ratios_CRC01_chr16.dummy.pkl.zip"
crc11_path = r"H:\Study\university\Computational-Biology\Year " \
             r"3\Projects\proj_scwgbs\resource\cpg_format\filtered_by_bl_and_cpgi\CRC11" \
             r"\all_cpg_ratios_CRC11_chr16.dummy.pkl.zip"
crc13_path = r"H:\Study\university\Computational-Biology\Year " \
             r"3\Projects\proj_scwgbs\resource\cpg_format\filtered_by_bl_and_cpgi\CRC13" \
             r"\all_cpg_ratios_CRC13_chr16.dummy.pkl.zip"

crc02_path = r"H:\Study\university\Computational-Biology\Year " \
             r"3\Projects\proj_scwgbs\resource\cpg_format\filtered_by_bl_and_cpgi\CRC02" \
             r"\all_cpg_ratios_CRC02_chr16.dummy.pkl.zip"
crc04_path = r"H:\Study\university\Computational-Biology\Year " \
             r"3\Projects\proj_scwgbs\resource\cpg_format\filtered_by_bl_and_cpgi\CRC04" \
             r"\all_cpg_ratios_CRC04_chr16.dummy.pkl.zip"
crc09_path = r"H:\Study\university\Computational-Biology\Year " \
             r"3\Projects\proj_scwgbs\resource\cpg_format\filtered_by_bl_and_cpgi\CRC09" \
             r"\all_cpg_ratios_CRC09_chr16.dummy.pkl.zip"

crc10_path = r"H:\Study\university\Computational-Biology\Year " \
             r"3\Projects\proj_scwgbs\resource\cpg_format\filtered_by_bl_and_cpgi\CRC10" \
             r"\all_cpg_ratios_CRC10_chr16.dummy.pkl.zip"
crc12_path = r"H:\Study\university\Computational-Biology\Year " \
             r"3\Projects\proj_scwgbs\resource\cpg_format\filtered_by_bl_and_cpgi\CRC12" \
             r"\all_cpg_ratios_CRC12_chr16.dummy.pkl.zip"
crc14_path = r"H:\Study\university\Computational-Biology\Year " \
             r"3\Projects\proj_scwgbs\resource\cpg_format\filtered_by_bl_and_cpgi\CRC14" \
             r"\all_cpg_ratios_CRC14_chr16.dummy.pkl.zip"

crc15_path = r"H:\Study\university\Computational-Biology\Year " \
             r"3\Projects\proj_scwgbs\resource\cpg_format\filtered_by_bl_and_cpgi\CRC15" \
             r"\all_cpg_ratios_CRC15_chr16.dummy.pkl.zip"


valid_path = r"H:\Study\university\Computational-Biology\Year 3\Projects\proj_scwgbs\covariance\valid_cpg.pkl"

if __name__ == '__main__':
    valid_data = pd.read_pickle(valid_path)
    valid_data = valid_data[valid_data["chromosome"] == "16"]
    valid_data["small_seq"] = valid_data["sequence"].str[73:77]

    cpg1 = valid_data[valid_data["sequence"].str.count("CG") == 1]
    cpg1["context"] = "other"

    cpg1.loc[cpg1["small_seq"].str.contains("[AT]CG[AT]", regex=True), "context"] = "WCGW"
    cpg1.loc[cpg1["small_seq"].str.contains("[CG]CG[CG]", regex=True), "context"] = "SCGS"

    only_needed = cpg1[["small_seq", "sequence", "context"]]
    only_needed = only_needed.transpose()
    only_needed.to_csv("info.csv")

    # crc01 = pd.read_pickle(crc01_path)
    # good = crc01[cpg1["location"]]
    # good.to_csv("crc01.csv")
    #
    # crc11 = pd.read_pickle(crc11_path)
    # good = crc11[cpg1["location"]]
    # good.to_csv("crc11.csv")
    #
    # crc13 = pd.read_pickle(crc13_path)
    # good = crc13[cpg1["location"]]
    # good.to_csv("crc13.csv")

    # rows = good.index.values
    # columns = list(good.columns.values)
    # data = good.values
    # data_added = np.vstack((data, cpg1["small_seq"]))
    # data_added = np.vstack((data_added, cpg1["context"]))
    # df = pd.DataFrame(data=data_added, index=columns + ["small_seq", "context"], columns=columns)

    crc02 = pd.read_pickle(crc02_path)
    good = crc02[cpg1["location"]]
    good.to_csv("crc02.csv")

    crc04 = pd.read_pickle(crc04_path)
    good = crc04[cpg1["location"]]
    good.to_csv("crc04.csv")

    crc09 = pd.read_pickle(crc09_path)
    good = crc09[cpg1["location"]]
    good.to_csv("crc09.csv")

    crc10 = pd.read_pickle(crc10_path)
    good = crc10[cpg1["location"]]
    good.to_csv("crc10.csv")

    crc12 = pd.read_pickle(crc12_path)
    good = crc12[cpg1["location"]]
    good.to_csv("crc12.csv")

    crc14 = pd.read_pickle(crc14_path)
    good = crc14[cpg1["location"]]
    good.to_csv("crc14.csv")

    crc15 = pd.read_pickle(crc15_path)
    good = crc15[cpg1["location"]]
    good.to_csv("crc15.csv")
