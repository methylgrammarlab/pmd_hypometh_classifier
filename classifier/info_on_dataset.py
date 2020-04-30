import matplotlib.pyplot as plt

from classifier.baseline_classifier import format_args, get_data

plt.style.use('seaborn-deep')


def split_wcgw_scgw_based_on_label(df, name):
    strong_rows = df[df["seq4"].str.contains("[CG]CG[CG]", regex=True)]
    weak_rows = df[df["seq4"].str.contains("[AT]CG[AT]", regex=True)]

    _ = plt.hist([strong_rows["label"], weak_rows["label"]], label=["SCGS", "WCGW"], normed=True)
    plt.xticks([0, 1], ["partial", "total"])
    plt.title("Histogram %s normed" % name)
    plt.legend(loc="upper right")
    plt.savefig("hist_%s_normed.png" % name)
    plt.close()

    _ = plt.hist([strong_rows["label"], weak_rows["label"]], label=["SCGS", "WCGW"])
    plt.xticks([0, 1], ["partial", "total"])
    plt.title("Histogram %s abs" % name)
    plt.legend(loc="upper right")
    plt.savefig("hist_%s_abs.png" % name)
    plt.close()


def main():
    args = format_args()
    train_data, test_data = get_data(args.data_path)
    split_wcgw_scgw_based_on_label(train_data, "CRC01")
    split_wcgw_scgw_based_on_label(test_data, "CRC13 and CRC11")


if __name__ == '__main__':
    main()
