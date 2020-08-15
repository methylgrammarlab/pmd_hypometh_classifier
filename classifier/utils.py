"""
Utils functions for the NN handling
Code adopted from  https://github.com/ohlerlab/DeepRiPe
"""

import glob
import os
import pickle
import re

import keras.backend as K
import numpy as np
import tensorflow as tf
from sklearn.metrics import roc_curve, auc, recall_score
from sklearn.model_selection import train_test_split, StratifiedKFold
from tensorflow.python.keras.models import load_model

SMALL_SEQ = "seq10"
BIG_SEQ = "sequence"

MOTIF_RE = re.compile("_*")
NO_KFOLD = 1


def get_train_test_data(path_to_data):
    """
    Get the test and train data from the pickle file, this function was created to allow other people to
    work with different data files
    :param path_to_data: The path for the data
    :return: the train and test dataset, each one is df
    """
    with open(path_to_data, "rb") as path_to_data_h:
        data = pickle.load(path_to_data_h)

    train_data, test_data = data["train"], data["test"]

    return train_data, test_data


def load_train_validate_test_data(path_to_data, input_len=150, only_test=False, validate_perc=0.2,
                                  kfold=5):
    """
    Load the train validate and test data and split it as requesetd using folds
    :param kfold: Number of folds to use
    :param validate_perc: The percentage of data used for validation
    :param only_test: Should we only extract the test data or all of it
    :param input_len: The input length
    :param path_to_data: path to file (consist of train, valid and test data)
    """
    x_train_list = []
    y_train_list = []
    x_validate_list = []
    y_validate_list = []

    if input_len == 150:
        seq_label = BIG_SEQ
    elif input_len == 10:
        seq_label = SMALL_SEQ
    else:
        raise NotImplemented("Unknown label to use")

    train_data, test_data = get_train_test_data(path_to_data)

    X_test_seq = np.array([seq_to_mat(seq) for seq in test_data[seq_label]])
    y_test = test_data["label"].values

    if only_test:
        return None, None, None, None, X_test_seq, y_test

    X_train_seq = np.array([seq_to_mat(seq) for seq in train_data[seq_label]])
    y_train = train_data["label"].values

    if kfold == NO_KFOLD:
        X_train_seq, X_valid_seq, y_train, y_valid = train_test_split(X_train_seq, y_train,
                                                                      test_size=validate_perc,
                                                                      random_state=42)
        x_train_list.append(X_train_seq)
        y_train_list.append(y_train)
        x_validate_list.append(X_valid_seq)
        y_validate_list.append(y_valid)

    else:
        kf = StratifiedKFold(n_splits=kfold, random_state=42, shuffle=True)
        for train_index, validation_index in kf.split(X_train_seq, y_train):
            X_train_fold, X_valid_fold = X_train_seq[train_index], X_train_seq[validation_index]
            y_train_fold, y_valid_fold = y_train[train_index], y_train[validation_index]

            x_train_list.append(X_train_fold)
            y_train_list.append(y_train_fold)
            x_validate_list.append(X_valid_fold)
            y_validate_list.append(y_valid_fold)

    return x_train_list, y_train_list, x_validate_list, y_validate_list, X_test_seq, y_test


########################
### custume metrics ####
########################
def accuracy(y_true, y_pred):
    y_true = tf.convert_to_tensor(y_true, np.float32)
    y_pred = tf.convert_to_tensor(y_pred, np.float32)

    diff = y_true - y_pred
    tptn = np.sum(diff == 0)
    return tptn / y_true.shape[0]


def precision(y_true, y_pred):
    y_true = tf.convert_to_tensor(y_true, np.float32)
    y_pred = tf.convert_to_tensor(y_pred, np.float32)

    TPs = K.sum(K.round(K.clip(y_true * y_pred, 0, 1)))
    predicted_positives = K.sum(K.round(K.clip(y_pred, 0, 1)))
    precision = TPs / (predicted_positives + K.epsilon())
    return precision


def precision_N(y_true, y_pred):
    y_true = tf.convert_to_tensor(y_true, np.float32)
    y_pred = tf.convert_to_tensor(y_pred, np.float32)
    y_true, y_pred = 1 - y_true, 1 - y_pred
    TNs = K.sum(K.round(K.clip(y_true * y_pred, 0, 1)))
    predicted_n = K.sum(K.round(K.clip(y_pred, 0, 1)))
    precision = TNs / (predicted_n + K.epsilon())
    return precision


def recall_TP(y_true, y_pred):
    y_true = tf.convert_to_tensor(y_true, np.float32)
    y_pred = tf.convert_to_tensor(y_pred, np.float32)
    TPs = K.sum(K.round(K.clip(y_true * y_pred, 0, 1)))
    possible_positives = K.sum(K.round(K.clip(y_true, 0, 1)))
    recall = TPs / (possible_positives + K.epsilon())
    return recall


def recall_TN(y_true, y_pred):
    y_true, y_pred = 1 - y_true, 1 - y_pred
    y_true = tf.convert_to_tensor(y_true, np.float32)
    y_pred = tf.convert_to_tensor(y_pred, np.float32)
    TNs = K.sum(K.round(K.clip(y_true * y_pred, 0, 1)))
    possible_positives = K.sum(K.round(K.clip(y_true, 0, 1)))
    recall = TNs / (possible_positives + K.epsilon())
    return recall


##########################################################
### Function to convert the sequence to one-hot encode ###
##########################################################

def seq_to_mat(seq):
    seq_len = len(seq)
    seq = seq.replace('A', '0')
    seq = seq.replace('a', '0')
    seq = seq.replace('C', '1')
    seq = seq.replace('c', '1')
    seq = seq.replace('G', '2')
    seq = seq.replace('g', '2')
    seq = seq.replace('T', '3')
    seq = seq.replace('t', '3')
    seq = seq.replace('U', '3')
    seq = seq.replace('u', '3')
    seq = seq.replace('N', '4')  # some cases have N in sequence
    seq = seq.replace('n', '4')
    seq_code = np.zeros((4, seq_len), dtype='float16')

    for i in range(seq_len):
        if int(seq[i]) != 4:
            seq_code[int(seq[i]), i] = 1

        else:
            seq_code[0:4, i] = np.tile(0.25, 4)

    return np.transpose(seq_code)


######################################################################################
######### function to find the top k and bottom K frequent 6MERs functions############
######################################################################################

def getkmer(X, y, pred, RBP_index, k):
    from rpy2.robjects.packages import importr
    base = importr('base')
    Biostrings = importr("Biostrings")

    multi_ind_high = [i[0] for i in sorted(enumerate(pred[:, RBP_index]), key=lambda x: x[1], reverse=True) if
                      y[i[0], RBP_index] == 1][0:k]
    multi_ind_low = [i[0] for i in sorted(enumerate(pred[:, RBP_index]), key=lambda x: x[1]) if
                     y[i[0], RBP_index] == 1][0:k]

    multi_fastaseq_low = vecs2dna(np.transpose(X[multi_ind_low], axes=(0, 2, 1)))
    multi_fastaseq_high = vecs2dna(np.transpose(X[multi_ind_high], axes=(0, 2, 1)))

    multi_fastaseq_high = base.unlist(multi_fastaseq_high)
    multi_fastaseq_low = base.unlist(multi_fastaseq_low)
    kmer_freqs_low = base.rowSums(
        base.sapply(Biostrings.DNAStringSet(multi_fastaseq_low),
                    Biostrings.oligonucleotideFrequency, width=6, step=1))
    kmer_freqs_high = base.rowSums(
        base.sapply(Biostrings.DNAStringSet(multi_fastaseq_high),
                    Biostrings.oligonucleotideFrequency, width=6, step=1))

    return kmer_freqs_low, kmer_freqs_high


###############################################################################
######### function to convert one hot encoded sequence to sequence ############
###############################################################################

def vecs2dna(seq_vecs):
    """
    Convert a list of sequences vectors to a sequence
    :param seq_vecs: np.array of sequences as one hot encoded
    :return: A list of sequences
    """
    trained_seq = []
    for i in range(0, seq_vecs.shape[0], 10000):
        temp = seq_vecs[i:i + 10000].astype(str)
        temp[np.all(temp == np.array(["1.0", "0.0", "0.0", "0.0"]), axis=2)] = "A"
        temp[np.all(temp == np.array(["0.0", "1.0", "0.0", "0.0"]), axis=2)] = "C"
        temp[np.all(temp == np.array(["0.0", "0.0", "1.0", "0.0"]), axis=2)] = "G"
        temp[np.all(temp == np.array(["0.0", "0.0", "0.0", "1.0"]), axis=2)] = "T"
        temp[np.all(temp == np.array(["0.25", "0.25", "0.25", "0.25"]), axis=2)] = "N"
        trained_seq.extend(["".join(i) for i in temp[:, :, 0]])

    return trained_seq


###########################################################
### Function to load kfold model and get scores from it ###
###########################################################
def load_models(models_folder):
    """
    Load the different models from the folder
    :param models_folder: path to the folder with models
    """
    models_paths = glob.glob(os.path.join(models_folder, "*.h5"))

    models = [load_model(model_path, custom_objects={'recall_TP': recall_TP, 'recall_TN': recall_TN}) for
              model_path in models_paths]
    return models


def get_scores(models, x_test, y_test):
    """
    Calculate and print the scores (accuracy and loss) for every model and the accuracy value
    :param models: A list of loaded models
    :param x_test: The dataset
    :param y_test: The real labels
    :return A list of accuracy scores and loss scores per fold
    """
    acc_per_fold = []
    loss_per_fold = []

    for model in models:
        score = model.evaluate(x_test, y_test)
        acc_per_fold.append(score[1] * 100)
        loss_per_fold.append(score[0])

    y_pred = predict(models, x_test)
    fpr_keras, tpr_keras, _ = roc_curve(y_test, y_pred)
    auc_keras = auc(fpr_keras, tpr_keras)
    recall = recall_score(y_test, np.round(y_pred))

    majority_vote_accuracy = np.sum((y_test - np.round(y_pred)) == 0) / y_test.shape[0]

    print('Average scores for all folds:')
    print(f'> Accuracy: {np.mean(acc_per_fold)} (+- {np.std(acc_per_fold)})')
    print(f'> Loss: {np.mean(loss_per_fold)}')

    print(f'>Truee accuracy using majority vote: {majority_vote_accuracy}')
    print(f'> AUC: {auc_keras}')
    print(f'> Recall: {recall}')


def predict(models, x):
    """
    Predict the output of x using the different models that were provided
    :param models: A list of loaded models
    :param x: The dataset to predict on
    :return the prediction for labels
    """
    y_pred = np.zeros(shape=(x.shape[0], len(models)))

    for i in range(len(models)):
        model_prediciton = models[i].predict(x)
        y_pred[:, i] = model_prediciton.reshape(model_prediciton.shape[0])

    return np.mean(y_pred, axis=1)
