import os

import numpy as np

np.random.seed(7)  # for reproducibility

import matplotlib.pyplot as plt
from concise.utils.plot import seqlogo_fig

from keras.models import load_model
from classifier.DeepRiPe.Scripts.IntegratedGradients import *
from classifier.DeepRiPe.Scripts.util_funcs import *
from classifier.DeepRiPe.Scripts.plotseqlogo import seqlogo_fig


#
# from tensorflow.python.keras.models import Model, load_model
# from classifier.utils import  load_data_merged
# from classifier.plotseqlogo import seqlogo_fig, seqlogo
# from classifier.IntegratedGradients import *


def theirs():
    # load data
    path_to_data = r"H:\Study\university\Computational-Biology\Year 3\Projects\proj_scwgbs\classifier\DeepRiPe\Data"

    path_to_datalow = os.path.join(path_to_data, "data_RBPsmed.h5")

    X_test_seq_low, X_test_region_low, y_test_RBP_low, y_test_name_low, y_train_low = load_data(
        path_to_datalow)

    # load models and obtain prediction and integrated_gradients

    path_to_model = r"H:\Study\university\Computational-Biology\Year " \
                    r"3\Projects\proj_scwgbs\classifier\DeepRiPe\Results\PARCLIP_models"
    path_to_modellow = os.path.join(path_to_model, "model_RBPsmed.h5")

    model_low = load_model(path_to_modellow, custom_objects={'precision': precision, 'recall': recall})
    pred_low = model_low.predict([X_test_seq_low, X_test_region_low])
    igres_low = integrated_gradients(model_low)

    # RBPnames for each model
    RBPnames_med = np.array(
        ['TARDBP', 'ELAVL2', 'ELAVL3', 'ELAVL4', 'RBM20', 'IGF2BP1', 'IGF2BP2', 'IGF2BP3', 'EWSR1', 'HNRNPD',
         'RBPMS', 'SRRM4', 'AGO2', 'NUDT21', 'FIP1L1', 'CAPRIN1', 'FMR1iso7', 'FXR2', 'AGO1', 'L1RE1',
         'ORF1'])

    # number of tasks for each model
    num_task_low = len(RBPnames_med)

    # one example from the low models, you can change the RBPname to see the attribution map of other RBPs in low models

    X_test_seq = X_test_seq_low
    X_test_region = X_test_region_low
    y_test_RBP = y_test_RBP_low

    RBPnames = RBPnames_med
    pred = pred_low
    igres = igres_low

    RBPname = "TARDBP"
    RBP_index = np.where(RBPnames == RBPname)[0][0]
    ind = []
    for i in sorted(enumerate(pred[:, RBP_index]), key=lambda x: x[1], reverse=True):
        if y_test_RBP[i[0], RBP_index] == 1 and pred[i[0], RBP_index] > 0.50:
            ind.append(i[0])

    # ind = ind[0:3]
    ex_seq = np.array(
        [igres.explain([X_test_seq[i], X_test_region[i]], outc=RBP_index, reference=False)[0] for i in ind])

    plt.close("all")
    seqlogo_fig(np.transpose(ex_seq[:, 25:125, :4], axes=(1, 2, 0)), vocab="RNA", figsize=(8, 3), ncol=1)
    plt.show()


def ours():
    x_train_seq, y_train, x_valid_seq, y_valid, x_test_seq, y_test = load_data_merged(
        r"dataset/classifier_data_ccpg1.pkl", 150, True)
    # model = load_model("test.h5",custom_objects={'precision': precision,'recall': recall })
    model = load_model(r"dataset/test.h5")
    pred = model.predict([x_test_seq])
    predict_label = 1
    ind = []

    for i in sorted(enumerate(pred), key=lambda x: x[1], reverse=True):
        if y_test[i[0]] == predict_label and pred[i[0]] > 0.50:
            ind.append(i[0])

    gradients = integrated_gradients(model)

    # Only first 10
    ind = ind[0:100]
    ex_seq = np.array([gradients.explain([x_test_seq[i]], reference=False)[0] for i in ind])


if __name__ == '__main__':
    theirs()
    # ours()
