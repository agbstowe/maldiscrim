################################################################################
####### Functional Neural Networks for MALDI-TOF spectrum classification #######
####### Creation date : 09-05-2025 #############################################
################################################################################
####### Usage : Data loading ###################################################
################################################################################

################################################################################
####### Import
################################################################################
import numpy as np
from typing import Optional
import os
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"
import tensorflow as tf
from sklearn.model_selection import train_test_split
from collections import Counter
import random
from itertools import accumulate

from .convolution import FunctionalConvolution
from .dense import FunctionalDense

def load_data(
        filepath: str,
        ratio: float = 0.8,
        n_labels: int = 1,
        batch_size: int = 32
) -> tuple:
    """
    Load data as tuple of TensorFlow datasets (split into training and test
    data).

    :param filepath: Path to csv file containing the data. The file should
        contain samples as rows and features as columns, where the last column
        corresponds to the label.
    :param ratio: Ratio of train and test data.
    :param n_labels: Number of label columns.
    :param batch_size: Size of batches.
    :return: TensorFlow dataset yielding labeled data of shape (24,), ().
    """
    data = tf.io.read_file(filepath)
    # data = tf.strings.split(data, sep='\r\n')[1:-1]
    data = tf.strings.split(data, sep='\n')[1:-1]
    data = tf.strings.split(data, sep=',').to_tensor()
    data = tf.strings.to_number(data)

    data, label = data[:, :-n_labels], data[:, -n_labels:]
    n_samples = data.shape[0]
    indices = list(range(n_samples))
    
    X_train, X_test, y_train, y_test = train_test_split(indices,
                                                        np.argmax(label, axis=1),
                                                        test_size=1-ratio,
                                                        stratify=np.argmax(label, axis=1),
                                                        random_state=42)
    datasets = []

    for idx in [X_train, X_test]:
        data_idx = tf.gather(data, idx)
        label_idx = tf.gather(label, idx)

        ds = tf.data.Dataset.from_tensor_slices((data_idx, label_idx))
        ds = ds.batch(batch_size)
        ds = ds.map(lambda x, y: (tf.expand_dims(x, -1), y))
        datasets.append(ds)

    train_data, test_data = datasets
    train_data = train_data.repeat(-1)

    return train_data, test_data


import random
from itertools import accumulate

def train_test(n_samples, train_size):
    """
    Sélectionne au hasard `train_size` indices dans chaque bloc défini par n_samples,
    puis renvoie la liste triée de tous ces indices.

    Parameters
    ----------
    n_samples : list of int
        Liste des tailles de chaque sous-population (ou bloc).
    train_size : int
        Nombre d'indices à prélever dans chaque bloc.

    Returns
    -------
    list of int
        Indices d'entraînement triés (0-based).
    """
    intervals = [0] + list(accumulate(n_samples))
    
    train_inds = []

    for i in range(len(n_samples)):
        start, end = intervals[i], intervals[i+1]
        block_inds = random.sample(range(start, end), train_size)
        train_inds.extend(block_inds)
    
    # return sorted(train_inds)
    return train_inds


def load_data_cv(
        filepath: str,
        train_size: int = 5,
        n_labels: int = 1,
        batch_size: int = 32
) -> tuple:
    """
    Load data as tuple of TensorFlow datasets (split into training and test
    data).

    :param filepath: Path to csv file containing the data. The file should
        contain samples as rows and features as columns, where the last column
        corresponds to the label.
    :param train_size: Number of samples of each category in train dataset.
    :param n_labels: Number of label columns.
    :param batch_size: Size of batches.
    :return: TensorFlow dataset yielding labeled data of shape (24,), ().
    """
    data = tf.io.read_file(filepath)
    # data = tf.strings.split(data, sep='\r\n')[1:-1]
    data = tf.strings.split(data, sep='\n')[1:-1]
    data = tf.strings.split(data, sep=',').to_tensor()
    data = tf.strings.to_number(data)

    data, label = data[:, :-n_labels], data[:, -n_labels:]
    n_samples = data.shape[0]
    indices = list(range(n_samples))
    
    n_samples = Counter(np.argmax(label, axis=1))

    X_train = train_test(n_samples=n_samples, train_size=train_size)
    X_test = np.setdiff1d(indices, X_train)

    datasets = []

    for idx in [X_train, X_test]:
        data_idx = tf.gather(data, idx)
        label_idx = tf.gather(label, idx)

        ds = tf.data.Dataset.from_tensor_slices((data_idx, label_idx))
        ds = ds.batch(batch_size)
        ds = ds.map(lambda x, y: (tf.expand_dims(x, -1), y))
        datasets.append(ds)

    train_data, test_data = datasets
    train_data = train_data.repeat(-1)

    return train_data, test_data

def setup_model(
        input_shape: tuple,
        filter_options: list,
        layer_options: list,
        loss: str = 'categorical_crossentropy',
        metrics: Optional[list] = None
) -> tf.keras.Model:
    """
    Setup model depending on the inputs shape and prediction aim. After a
    channel-wise normalization, FunctionalConvolution and FunctionalDense
    layers are added (according to the specifications).

    :param input_shape: Shape of input data (excluding batch dimension), e.g.
        - (time, n_channels) for multivariate time series
        - (space1, space2, n_channels) for images
        - (space1, space2, time, n_channels) for videos
    :param filter_options: List of dictionaries specifying the options used to
        create FunctionalConvolution layers. The number of convolutional layers
        equals the length of the list.
    :param layer_options: List of dictionaries specifying the options used to
        create FunctionalDense layers. The number of dense layers equals the
        length of the list.
    :param loss: Loss function to be used during training.
    :param metrics: List of metrics to be used during training.
    :return: Compiled model
    """
    if metrics is None:
        metrics = ['accuracy']

    inputs = tf.keras.layers.Input(shape=input_shape)

    # norm_axes = list(range(len(input_shape) - 1))
    # layer = tf.keras.layers.LayerNormalization(
    #     axis=-1,
    #     center=False,
    #     scale=False,
    #     epsilon=1e-10,
    #     name='Normalization'
    # )(inputs)

    a = 0
    for i, filter_option in enumerate(filter_options):
      if a==0:
        layer = FunctionalConvolution(
            **filter_option,
            name=f'FunctionalConvolution_{i}'
        )(inputs)
        a += 1
      else:
        layer = FunctionalConvolution(
            **filter_option,
            name=f'FunctionalConvolution_{i}'
        )(layer)

    for i, layer_option in enumerate(layer_options):
        layer = FunctionalDense(
            **layer_option,
            name=f'FunctionalDense_{i}'
        )(layer)

    outputs = layer

    model = tf.keras.Model(inputs=inputs, outputs=outputs, name="Example_FNN")
    model.compile(loss=loss, optimizer='adam', metrics=metrics)

    return model
