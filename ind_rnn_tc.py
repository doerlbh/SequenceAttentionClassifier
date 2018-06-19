from models.modules.indRNN import IndRNNCell
from models.attention import attention
from tensorflow.contrib.rnn import BasicLSTMCell
from tensorflow.python.ops.rnn import bidirectional_dynamic_rnn as bi_rnn
import time
from utils.prepare_data import *

# Hyperparameter
MAX_DOCUMENT_LENGTH = 100
EMBEDDING_SIZE = 128
HIDDEN_SIZE = 64
ATTENTION_SIZE = 64
lr = 1e-3
BATCH_SIZE = 256
KEEP_PROB = 0.5
LAMBDA = 1e-3
MAX_LABEL = 2
epochs = 20

vocab_size = 5
MAX_SEQ_LENGTH = MAX_DOCUMENT_LENGTH

vocab = ["N","A","C","T","G"]

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

# split dataset to test and dev
x_test, x_dev, y_test, y_dev, dev_size, test_size = \
    split_dataset(x_test, y_test, 0.1)
print("Validation size: ", dev_size)

graph = tf.Graph()
with graph.as_default():

    batch_x = tf.placeholder(tf.int32, [None, MAX_DOCUMENT_LENGTH])
    batch_y = tf.placeholder(tf.float32, [None, MAX_LABEL])
    keep_prob = tf.placeholder(tf.float32)

    embeddings_var = tf.Variable(tf.random_uniform([vocab_size, EMBEDDING_SIZE], -1.0, 1.0), trainable=True)
    batch_embedded = tf.nn.embedding_lookup(embeddings_var, batch_x)
    print(batch_embedded.shape)  # (?, 256, 100)

    cell = IndRNNCell(HIDDEN_SIZE)
    rnn_outputs, _ = tf.nn.dynamic_rnn(cell, batch_embedded, dtype=tf.float32)

    # Attention
    attention_output, alphas = attention(rnn_outputs, ATTENTION_SIZE, return_alphas=True)
    drop = tf.nn.dropout(attention_output, keep_prob)
    shape = drop.get_shape()

    # Fully connected layer（dense layer)
    W = tf.Variable(tf.truncated_normal([shape[1].value, MAX_LABEL], stddev=0.1))
    b = tf.Variable(tf.constant(0., shape=[MAX_LABEL]))
    y_hat = tf.nn.xw_plus_b(drop, W, b)

    loss = tf.reduce_mean(tf.nn.sigmoid_cross_entropy_with_logits(logits=y_hat, labels=batch_y))
    optimizer = tf.train.AdamOptimizer(learning_rate=lr).minimize(loss)

    # Accuracy metric
    prediction = tf.argmax(tf.nn.softmax(y_hat), 1)
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
        print("Validation accuracy: ", sess.run([accuracy, loss], feed_dict={
            batch_x: x_dev,
            batch_y: y_dev,
            keep_prob: 1.0
        }))

    print("Training finished, time consumed : ", time.time() - start, " s")
    print("start predicting:  \n")
    test_accuracy = sess.run([accuracy], feed_dict={batch_x: x_test, batch_y: y_test, keep_prob: 1})
    print("Test accuracy : %f %%" % (test_accuracy[0] * 100))

    save_path = saver.save(sess, "./ind_rnn_tc.ckpt")
    print("Model saved in path: %s" % save_path)





