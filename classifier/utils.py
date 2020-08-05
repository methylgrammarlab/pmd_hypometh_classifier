"""
Utils functions for the NN handling
Code adopted from  https://github.com/ohlerlab/DeepRiPe
"""

import pickle
import re

import keras.backend as K
import numpy as np
import tensorflow as tf
from sklearn.model_selection import train_test_split, StratifiedKFold

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


def load_train_validate_test_data(path_to_data, input_len=150, only_test=False, validate_perc=0.2, kfold=5):
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


###############################################################################
######### function to convert one hot encoded sequence to sequence ############
###############################################################################

def vecs2dna(seq_vecs):
    """
    Convert a list of sequences vectors to a sequence
    :param seq_vecs: np.array of sequences as one hot encoded
    :return: A list of sequences
    """
    seq_vecs = seq_vecs.astype(str)
    seq_vecs[np.all(seq_vecs == np.array(["1.0", "0.0", "0.0", "0.0"]), axis=2)] = "A"
    seq_vecs[np.all(seq_vecs == np.array(["0.0", "1.0", "0.0", "0.0"]), axis=2)] = "C"
    seq_vecs[np.all(seq_vecs == np.array(["0.0", "0.0", "1.0", "0.0"]), axis=2)] = "G"
    seq_vecs[np.all(seq_vecs == np.array(["0.0", "0.0", "0.0", "1.0"]), axis=2)] = "T"
    seq_vecs[np.all(seq_vecs == np.array(["0.25", "0.25", "0.25", "0.25"]), axis=2)] = "N"
    sequences = ["".join(i) for i in seq_vecs[:, :, 0]]
    return sequences
