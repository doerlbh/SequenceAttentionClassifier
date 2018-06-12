import time

from models.modules.multihead import *
from utils.prepare_data import *

import pandas as pd

# Hyperparameter
MAX_SEQ_LENGTH = 100
EMBEDDING_SIZE = 10
HIDDEN_SIZE = 512
ATTENTION_SIZE = 64
lr = 1e-3
BATCH_SIZE = 256
KEEP_PROB = 0.5
LAMBDA = 0.0001

MAX_LABEL = 2

vocab_size = 5

epochs = 5

# load data
x_train, y_train = load_data("./data/train-BRAF.csv", sample_ratio=1)
x_test, y_test = load_data("./data/test-BRAF.csv", sample_ratio=1)

# data preprocessing
base_number = {'N':0, 'A':1, 'C':2, 'T':3, 'G':4}
x_train_l = np.ndarray((len(x_train),MAX_SEQ_LENGTH))

for t in np.arange(len(x_train)):
    line = list(x_train[t])[1:MAX_SEQ_LENGTH+1]
    for k in np.arange(MAX_SEQ_LENGTH):
        x_train_l[t,k] = base_number[line[k]]
    
x_test_l = np.ndarray((len(x_test),MAX_SEQ_LENGTH))

for t in np.arange(len(x_test)):    
    line = list(x_test[t])[1:MAX_SEQ_LENGTH+1]
    for k in np.arange(MAX_SEQ_LENGTH):
        x_test_l[t,k] = base_number[line[k]]
    
print(x_train_l.shape)
print(x_test_l.shape)

# x_train = x_train_l[0:100000,:]
# x_test = x_test_l[0:10000,:]
x_train = x_train_l
x_test = x_test_l

# y_train = y_train[0:100000,:]
# y_test = y_test[0:10000,:]

# split dataset to test and dev
x_test, x_dev, y_test, y_dev, dev_size, test_size = \
    split_dataset(x_test, y_test, 0.1)
print("Validation size: ", dev_size)

graph = tf.Graph()
with graph.as_default():

    batch_x = tf.placeholder(tf.int32, [None, MAX_SEQ_LENGTH])
    batch_y = tf.placeholder(tf.float32, [None, MAX_LABEL])
    keep_prob = tf.placeholder(tf.float32)

    embeddings_var = tf.Variable(tf.random_uniform([vocab_size, EMBEDDING_SIZE], -1.0, 1.0), trainable=True)
    batch_embedded = tf.nn.embedding_lookup(embeddings_var, batch_x)
    # multihead attention
    outputs = multihead_attention(queries=batch_embedded, keys=batch_embedded)
    # FFN(x) = LN(x + point-wisely NN(x))
    outputs = feedforward(outputs, [HIDDEN_SIZE, EMBEDDING_SIZE])
    print(outputs.shape)
    outputs = tf.reshape(outputs, [-1, MAX_SEQ_LENGTH * EMBEDDING_SIZE])
    logits = tf.layers.dense(outputs, units=MAX_LABEL)
    loss = tf.reduce_mean(tf.nn.sigmoid_cross_entropy_with_logits(logits=logits, labels=batch_y))
    optimizer = tf.train.AdamOptimizer(learning_rate=lr).minimize(loss)

    # Accuracy metric
    prediction = tf.argmax(tf.nn.softmax(logits), 1)
    accuracy = tf.reduce_mean(tf.cast(tf.equal(prediction, tf.argmax(batch_y, 1)), tf.float32))

with tf.Session(graph=graph) as sess:
    sess.run(tf.global_variables_initializer())
    print("Initialized! ")

    print("Start trainning")
    start = time.time()
    for e in range(epochs):

        epoch_start = time.time()
        print("Epoch %d start !" % (e + 1))
        for x_batch, y_batch in fill_feed_dict(x_train, y_train, BATCH_SIZE):
            fd = {batch_x: x_batch, batch_y: y_batch, keep_prob: KEEP_PROB}
            l, _, acc = sess.run([loss, optimizer, accuracy], feed_dict=fd)

        epoch_finish = time.time()
        print("Validation accuracy and loss: ", sess.run([accuracy, loss], feed_dict={
            batch_x: x_dev,
            batch_y: y_dev,
            keep_prob: 1.0
        }))
        print("epoch time:", epoch_finish - epoch_start , " s")

    print("Training finished, time consumed : ", time.time() - start, " s")
    print("start predicting:  \n")
    test_accuracy = sess.run([accuracy], feed_dict={batch_x: x_test, batch_y: y_test, keep_prob: 1})
    print("Test accuracy : %f %%" % (test_accuracy[0] * 100))




