import argparse
import os
import pickle
import sys

import numpy as np

np.random.seed(7)  # for reproducibility

from sklearn.model_selection import train_test_split

# tf.python.control_flow_ops = tf


from tensorflow.python.keras.models import Model, load_model
from tensorflow.python.keras.layers import Input
from tensorflow.python.keras.layers import Dense, Dropout
from tensorflow.python.keras.layers.convolutional import Conv1D
from tensorflow.python.keras.layers.pooling import MaxPooling1D
from tensorflow.python.keras.optimizers import Adam
from tensorflow.python.keras.callbacks import ModelCheckpoint, EarlyStopping
import tensorflow.python.keras.backend as K

import matplotlib as mpl

mpl.use('Agg')
# from keras.utils.layer_utils import print_layer_shapes

path_to_scripts = "./Scripts/"
sys.path.append(path_to_scripts)
from classifier.DeepRiPe.Scripts.util_funcs import precision, recall, seq_to_mat

SEQ_CONSTANT = "seq10"


def load_data_merged(path_to_data):
    """
    load the data
    :param path_to_data: path to file (consist of train, valid and test data)
    """
    with open(path_to_data, "rb") as path_to_data_h:
        data = pickle.load(path_to_data_h)

    train_data, test_data = data["train"], data["test"]

    X_train_seq = np.array([seq_to_mat(seq) for seq in train_data[SEQ_CONSTANT]])
    y_train = train_data["label"].values

    X_test_seq = np.array([seq_to_mat(seq) for seq in test_data[SEQ_CONSTANT]])
    y_test = test_data["label"].values

    # TODO: in the future need to see how to do it with cross-validation
    X_train_seq, X_valid_seq, y_train, y_valid = train_test_split(X_train_seq, y_train, test_size=0.2,
                                                                  random_state=42)

    return X_train_seq, y_train, X_valid_seq, y_valid, X_test_seq, y_test


def create_seq_model(input_len):
    """
    Create a sequence model
    :param input_len: path to file (consist of train, valid and test data)
    """
    K.clear_session()
    # tf.random.set_seed(5005)
    dim = 4  # TODO: this might need to be 1

    filters_num = [90, 100]
    filter_len = [7, 7]
    pool_len = [4, 10]
    strides = [2, 5]

    input_node = Input(shape=(input_len, dim), name="input")
    conv1 = Conv1D(filters=filters_num[0], kernel_size=filter_len[0], padding='valid', activation="relu",
                   name="conv1")(input_node)
    pool1 = MaxPooling1D(pool_size=pool_len[0], strides=strides[0], name="left_pool1")(conv1)
    drop1 = Dropout(0.25, name="left_drop1")(pool1)

    # conv_merged = Conv1D(filters=100, kernel_size=5, padding='valid', activation="relu", name="conv_merged")(
    #     drop1)
    # merged_pool = MaxPooling1D(pool_size=10, strides=5)(conv_merged)
    # merged_drop = Dropout(0.25)(merged_pool)
    # merged_flat = Flatten()(merged_drop)

    hidden1 = Dense(250, activation='relu', name="hidden1")(drop1)
    output = Dense(1, activation='sigmoid', name="output")(hidden1)
    model = Model(inputs=[input_node], outputs=output)
    print(model.summary())
    return model


def train_diff_model(data_path, res_path, model_name, input_len,
                     num_epoch, batchsize, model_path="./weights.hdf5"):
    """
    Training the model
    :param data_path: path to file (consist of train, valid and test data)
    :param res_path:
    :param model_name:
    :param input_len:
    :param num_epoch:
    :param batchsize:
    :param model_path:
    :return:
    """
    filter_lengths = [4, 5]
    print('creating model')
    model = create_seq_model(input_len)
    print('compiling model')
    adam = Adam(lr=0.001, beta_1=0.9, beta_2=0.999, epsilon=1e-08, decay=1e-6)
    model.compile(loss='binary_crossentropy', optimizer=adam, metrics=['accuracy', precision, recall])
    checkpointer = ModelCheckpoint(filepath=model_path, verbose=1, save_best_only=True)
    earlystopper = EarlyStopping(monitor='val_loss', patience=5, verbose=1)
    # tb=TensorBoard(log_dir='./Output/logs', histogram_freq=0, write_graph=True, write_images=False)

    print('loading data')
    x_train_seq, y_train, x_valid_seq, y_valid, x_test_seq, y_test = load_data_merged(data_path)

    print('fitting the model')
    history = model.fit(x_train_seq, y_train, epochs=num_epoch, batch_size=batchsize,
                        validation_data=(x_valid_seq, y_valid), verbose=2,
                        callbacks=[checkpointer, earlystopper, ])  # tb])

    print('saving the model')
    model.save(os.path.join(res_path, model_name + ".h5"))

    print('testing the model')
    score = model.evaluate(x_test_seq, y_test)

    print(model.metrics_names)
    for i in range(len(model.metrics_names)):
        print(str(model.metrics_names[i]) + ": " + str(score[i]))

    print("{}: {:.2f}".format(model.metrics_names[0], score[0]))
    print("{}: {:.2f}".format(model.metrics_names[1], score[1]))
    print("{}: {:.2f}".format(model.metrics_names[2], score[2]))


################################################################################
# Testing the model
#
# Input: path to file (consist of train, valid and test data)
#
################################################################################

def test_model(output_path, data_path, res_path, model_name):
    print('test the model and plot the curve')
    model = load_model(os.path.join(res_path, model_name + ".h5"),
                       custom_objects={'precision': precision, 'recall': recall})

    with open(data_path, "rb") as data_path_h:
        data = pickle.load(data_path_h)

    test_data = data["test"]
    X_test_seq = test_data[SEQ_CONSTANT].values
    y_test = test_data["label"].values

    print('predicting on test data')
    y_pred = model.predict(X_test_seq, verbose=1)
    model.evaluate(X_test_seq, y_test)

    print("saving the prediction to " + output_path)
    with open(output_path, "wb") as output:
        pickle.dump(y_pred, output)


def format_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-il', dest='input_len', default=None, type=int, help='input length')
    parser.add_argument('-ne', dest='num_epoch', default=None, type=int, help='number of epochs')
    parser.add_argument('-bs', dest='batchsize', default=None, type=int, help='Batch size')
    parser.add_argument('-dp', dest='data_path', default=None, type=str, help='path to the data')
    parser.add_argument('-op', dest='output_path', default=None, type=str, help='path to the output')
    parser.add_argument('-mp', dest='model_path', default=None, type=str, help='path to the model')
    parser.add_argument('-pp', dest='prediction_path', default=None, type=str, help='path to the prediction')
    parser.add_argument('-name', dest='model_name', default=None, type=str, help='name of the model')
    parser.add_argument('-t', dest='test', default=False, type=bool, help='test the model')
    return parser.parse_args()


def main():
    args = format_args()

    # train_diff_model(data_path=args.data_path, res_path=args.output_path, model_name=args.model_name,
    #                  input_len=args.input_len, num_epoch=args.num_epoch, batchsize=args.batchsize,
    #                  model_path=args.model_path)
    train_diff_model(
        data_path=r"H:\Study\university\Computational-Biology\Year 3\Projects\proj_scwgbs\classifier\classifier_data_ccpg1.pkl",
        res_path=".", model_name="test",
        input_len=10, num_epoch=100, batchsize=32,
        model_path="model")

    if args.test:
        print("testing the model and plot the curves")
        test_model(output_path=args.prediction_path, data_path=args.data_path, res_path=args.output_path,
                   model_name=args.model_name)


if __name__ == '__main__':
    main()
