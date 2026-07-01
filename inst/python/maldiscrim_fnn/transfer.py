################################################################################
####### Functional Neural Networks for MALDI-TOF spectrum classification #######
####### Creation date : 30-06-2026 #############################################
################################################################################
####### Usage : Transfer learning ###############################################
################################################################################

################################################################################
####### Import
################################################################################
from typing import Optional

import tensorflow as tf

from .convolution import FunctionalConvolution
from .dense import FunctionalDense


def freeze_layers(model: tf.keras.Model, n_freeze_layers: int) -> tf.keras.Model:
    """
    Freeze the first `n_freeze_layers` FunctionalConvolution/FunctionalDense
    layers of a loaded model (in network order), leaving every other layer
    trainable. The Input layer is not counted.

    :param model: A loaded Keras model (e.g. via .loadFNNModel()).
    :param n_freeze_layers: Number of leading functional layers to freeze.
        0 freezes nothing (full fine-tuning); a value greater than or equal
        to the number of functional layers in the model freezes them all.
    :return: The same model object, mutated in place (layers marked
        trainable/non-trainable) and returned for convenience.
    """
    functional_layers = [
        layer for layer in model.layers
        if isinstance(layer, (FunctionalConvolution, FunctionalDense))
    ]

    for i, layer in enumerate(functional_layers):
        layer.trainable = i >= n_freeze_layers

    return model


def setup_transfer_model(
        base_model: tf.keras.Model,
        n_freeze_layers: int,
        new_layer_options: list,
        loss: str = 'categorical_crossentropy',
        metrics: Optional[list] = None
) -> tf.keras.Model:
    """
    Build a transfer-learning model from a previously trained FNN.

    The base model is truncated just before its final FunctionalDense
    layer(s): the output of the last FunctionalConvolution layer (i.e. the
    last layer that is not a FunctionalDense) is reused as input to one or
    more freshly initialized FunctionalDense layers, specified by
    `new_layer_options` (same format as `layer_options` in `setup_model`).
    This allows the new head to target a different number of classes than
    the source model.

    Layers carried over from the base model are frozen/unfrozen according to
    `n_freeze_layers` (see `freeze_layers`); the newly created head layers
    are always trainable.

    :param base_model: A loaded Keras FNN model (e.g. via .loadFNNModel()),
        built with `setup_model` from data_model_loading.py.
    :param n_freeze_layers: Number of leading functional layers (in network
        order) of `base_model` to freeze. Use the number of
        FunctionalConvolution layers to freeze all of them
        (feature-extraction), one fewer to leave the last convolutional
        layer trainable (partial fine-tuning), or 0 to fine-tune the full
        carried-over backbone.
    :param new_layer_options: List of dictionaries specifying the options
        used to create the new FunctionalDense head layer(s). Same format
        as `layer_options` in `setup_model` (each dict accepts n_neurons,
        basis_options, activation, pooling).
    :param loss: Loss function to be used during training.
    :param metrics: List of metrics to be used during training.
    :return: A new, compiled tf.keras.Model with the transferred backbone
        and the new head, ready to be fit.
    """
    if metrics is None:
        metrics = ['accuracy']

    # Identify the last layer that is NOT a FunctionalDense head layer.
    # Everything up to and including that layer is reused as the backbone;
    # the original FunctionalDense layer(s) are discarded and replaced.
    backbone_output = None
    cut_index = None

    for i, layer in enumerate(base_model.layers):
        if not isinstance(layer, FunctionalDense):
            cut_index = i

    if cut_index is None:
        raise ValueError(
            "No non-FunctionalDense layer found in base_model; cannot "
            "determine where to cut the backbone for transfer learning."
        )

    backbone_output = base_model.layers[cut_index].output

    # Build a backbone-only model so we can freeze/unfreeze its layers
    # independently of the (discarded) original head, then reuse its
    # output tensor as the entry point of the new functional graph.
    backbone = tf.keras.Model(
        inputs=base_model.input,
        outputs=backbone_output,
        name="FNN_backbone"
    )
    backbone = freeze_layers(backbone, n_freeze_layers)

    layer = backbone.output
    for i, layer_option in enumerate(new_layer_options):
        layer = FunctionalDense(
            **layer_option,
            name=f'FunctionalDense_transfer_{i}'
        )(layer)

    outputs = layer

    model = tf.keras.Model(
        inputs=backbone.input,
        outputs=outputs,
        name="Transfer_FNN"
    )
    model.compile(loss=loss, optimizer='adam', metrics=metrics)

    return model
