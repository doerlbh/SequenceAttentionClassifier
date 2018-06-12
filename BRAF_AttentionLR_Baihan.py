
# coding: utf-8

# ### Notebook to configure model

# In[2]:


import time

import numpy as np
from models.modules.multihead import *
from utils.prepare_data import *

import pandas as pd


import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_context(context="talk")

import torch
import torch.nn as nn
import tensorflow as tf
import torch.nn.functional as F
from torchvision import datasets
import torchvision.transforms as transforms
from torch.autograd import Variable
from torch.utils import data

import math
import copy

# from keras.utils import np_utils


# In[3]:


# Hyperparameter

QUERY_SIZE = 64
EMBEDDING_SIZE = 128

# HIDDEN_SIZE = 512
# ATTENTION_SIZE = 64

LEARNING_RATE = 1e-3
BATCH_SIZE = 50

VOCAB_SIZE = 5
NUM_CLASSES = 2

# Data-specific

READ_LENGTH = 100

GENOME_START = 140719327
GENOME_END = 140924764

GENOME_LENGTH = GENOME_END - GENOME_START + 1
CONTEXT_SIZE = GENOME_LENGTH


# In[28]:


class TensorizedReadDataset(torch.utils.data.DataLoader):
    'Characterizes a Tensorized dataset for genome reads in PyTorch'
    
    def __init__(self, reads, ref_locs, labels, read_length=100, genome_start=0, genome_end=0):
#         super(TensorizedReadDataset, self).__init__()
        
        self.read_length = read_length
        self.labels = labels
        self.reads = reads
        self.ref_locs = ref_locs
        self.genome_start = genome_start
        self.genome_end = genome_end

    def __len__(self):
        return len(self.reads)

    def __getitem__(self, index):
        
        vals = list(self.reads[index])[0:self.read_length]        
        locs = list(np.arange(self.ref_locs[index]-self.genome_start,self.ref_locs[index]+self.read_length-self.genome_start))

#         print(len(vals))
#         print(len(locs))
        
        vals2idx = {'N': 0, 'A': 1, 'C': 2, 'T': 3, 'G': 4}
#         read = torch.LongTensor(np.array([vals2idx[val]+loc*len(vals2idx) for val, loc in zip(vals, locs)], dtype=int), requires_grad=False)

        read = torch.autograd.Variable(torch.LongTensor(np.array([vals2idx[val]+loc*len(vals2idx) for val, loc in zip(vals, locs)], dtype=int)), requires_grad=False)
        
        X = read
        Y = self.labels[index,:]

        return X, Y
    


# In[10]:


def attention(query, key, value, dropout=None):
    
    d_k = query.size(-1)
    scores = torch.matmul(query, key.transpose(-2, -1))              / math.sqrt(d_k)
    p_attn = F.softmax(scores, dim = -1)
    if dropout is not None:
        p_attn = dropout(p_attn)
    return torch.matmul(p_attn, value), p_attn


# In[11]:


class NGramDenseEmbedding(nn.Module):

    def __init__(self, vocab_size, embedding_dim, context_size):
        super(NGramDenseEmbedding, self).__init__()
        self.embeddings = nn.Embedding(vocab_size, embedding_dim)
        self.linear1 = nn.Linear(context_size * embedding_dim, 128)
        self.linear2 = nn.Linear(128, vocab_size)

    def forward(self, inputs):
        embeds = self.embeddings(inputs).view((1, -1))
        out = F.relu(self.linear1(embeds))
        out = self.linear2(out)
        log_probs = F.log_softmax(out, dim=1)
        return log_probs


# In[12]:


class AttentionLR(nn.Module):
    
    def __init__(self, input_size, num_classes, dropout=0.1):
        super(AttentionLR, self).__init__()
        
        self.dropout = nn.Dropout(p=dropout)
        self.KQ_attn = None
        self.KQV_attn = None
        self.linears = clones(nn.Linear(input_size, num_classes),1)
        
        self.K = NGramDenseEmbedding(VOCAB_SIZE, EMBEDDING_SIZE, CONTEXT_SIZE)
        self.V = NGramDenseEmbedding(VOCAB_SIZE, EMBEDDING_SIZE, CONTEXT_SIZE)
        
        self.linear = nn.Linear(input_size, num_classes)
        
    def forward(self, query_seq):
        
        Q_lookup = torch.tensor([word_to_ix[w] for w in query_seq], dtype=torch.long)
        
        K_lookup = self.K(Q_lookup)
        V_lookup = self.V(Q_lookup)       
        
        self.KQ_attn, self.KQV_attn = attention(Q_lookup, K_lookup, V_lookup, dropout=self.dropout)
                
        return F.log_softmax(self.linear(self.KQV_attn), dim=1)
    


# In[62]:


class SequenceAttentionClassifier(nn.Module):
    
    def __init__(self, genome_length, read_length=100, vocab_size=5, query_size=64, embedding_size=128, num_classes=2):
        super(SequenceAttentionClassifier, self).__init__()
        self.genome_length = genome_length
        self.read_length = read_length
        self.vocab_size = vocab_size
        self.query_size = query_size
        self.embedding_size = embedding_size
        self.num_classes = num_classes
        self.K = nn.Embedding(vocab_size*genome_length, embedding_size)
        self.V = nn.Embedding(vocab_size*genome_length, query_size)
        self.W = nn.Linear(query_size, num_classes)
        self.Q = nn.Linear(embedding_size, query_size)
        
    def forward(self, read):
        
        # 'read' here should be mapped to a flattened form where X_ij = 1 maps to i*vocab_size + j
        K_lookup = self.K(read) # Get the relevant keys
        V_lookup = self.V(read) # Get the relevant values

        # Get the attention weights
        logits = self.Q(K_lookup) / math.sqrt(self.embedding_size)
        probs = F.softmax(logits, dim = -1)
                
        # Calculate the covariates for the logistic regression
        X = torch.matmul(probs.transpose(1,2), V_lookup)

        # Right now we can just ignore the fact that we're doing a linear-transform.
        # In the future we'll add nonlinearities

        # Return the logits for the classifier
        return self.W(X)
    


# In[6]:


# load data
x_train, refs_train, y_train = load_ref_data("../data/ref-train-BRAF.csv", sample_ratio=1)
x_test, refs_test, y_test = load_ref_data("../data/ref-test-BRAF.csv", sample_ratio=1)

# split dataset to test and dev
x_train, x_softval, y_train, y_softval, refs_train, refs_softval, softval_size, train_size =     split_ref_dataset(x_train, y_train, refs_train, 0.01)
    
print("Soft Validation size: ", softval_size)
print("Training size: ", train_size)


# In[50]:


# Generators
train_dataset = TensorizedReadDataset(reads=x_train, 
                                      ref_locs=refs_train, 
                                      labels=y_train, 
                                      read_length=READ_LENGTH, 
                                      genome_start=GENOME_START, 
                                      genome_end=GENOME_END)

hardval_dataset = TensorizedReadDataset(reads=x_test, 
                                        ref_locs=refs_test, 
                                        labels=y_test, 
                                        read_length=READ_LENGTH, 
                                        genome_start=GENOME_START, 
                                        genome_end=GENOME_END)

softval_dataset = TensorizedReadDataset(reads=x_softval, 
                                        ref_locs=refs_softval, 
                                        labels=y_softval, 
                                        read_length=READ_LENGTH, 
                                        genome_start=GENOME_START, 
                                        genome_end=GENOME_END)

# Input pipeline
train_loader = torch.utils.data.DataLoader(dataset=train_dataset,
                                           batch_size=BATCH_SIZE,
                                           shuffle=True)

hardval_loader = torch.utils.data.DataLoader(dataset=hardval_dataset,
                                             batch_size=BATCH_SIZE,
                                             shuffle=True)

softval_loader = torch.utils.data.DataLoader(dataset=softval_dataset,
                                             batch_size=BATCH_SIZE,
                                             shuffle=True)


# In[63]:


# model = AttentionLR(MAX_SEQ_LENGTH, MAX_LABEL)
model = SequenceAttentionClassifier(genome_length=GENOME_LENGTH,
                                    read_length=READ_LENGTH,
                                    vocab_size=VOCAB_SIZE,
                                    query_size=QUERY_SIZE,
                                    embedding_size=EMBEDDING_SIZE,
                                    num_classes=NUM_CLASSES)

criterion = nn.CrossEntropyLoss()
optimizer = torch.optim.SGD(model.parameters(), lr=LEARNING_RATE)


# In[71]:


num_epochs = 1

# Training process

b = 0 # count batch
for epoch in range(num_epochs):
    for x_batch, y_batch in train_loader:
        
        optimizer.zero_grad()
        
        outputs = model(x_batch) 
        
#         for debugging the output size
        print(outputs.shape)
        print(y_batch.shape)
        print(type(outputs))
        print(type(y_batch))
                
        loss = criterion(outputs, y_batch)  
        loss.backward()
        optimizer.step()
        
        if (b + 1) % 10 == 0:
            print("Epoch {}, Batch {}, loss :{}".format(epoch + 1, b + 1, loss.data[0]))
        b = b + 1
            


# In[ ]:


# Hard Validation process

correct = 0
total = 0
for i, (x_batch, y_batch) in hardval_loader:
    outputs = model(x_batch)
    _, predicted = torch.max(outputs.data, 1)
    total += labels.size(0)
    correct += (predicted == labels).sum()

print('Hard validation Accuracy: {}%'.format(100 * correct / total))    


# Soft Validation process

correct = 0
total = 0
for i, (x_batch, y_batch) in softval_loader:
    outputs = model(x_batch)
    _, predicted = torch.max(outputs.data, 1)
    total += labels.size(0)
    correct += (predicted == labels).sum()

print('Soft validation Accuracy: {}%'.format(100 * correct / total))    

