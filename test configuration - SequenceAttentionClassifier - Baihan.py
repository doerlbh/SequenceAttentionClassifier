
# coding: utf-8

# ### Notebook to configure model

# In[1]:


import time
import math
import copy

import numpy as np
import pandas as pd

import matplotlib
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


# In[2]:


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

ref_names = ["class", "ref", "sequence"]


# In[4]:


def load_ref_data(file_name, sample_ratio= 1, n_class=2, names=ref_names):
    
    csv_file = pd.read_csv(file_name, names=ref_names)
    shuffle_csv = csv_file.sample(frac=sample_ratio).reset_index()
#     x = pd.Series(shuffle_csv["sequence"])
    x = list(shuffle_csv["sequence"])
#     ref = pd.Series(shuffle_csv["ref"])
    ref = list(shuffle_csv["ref"])
    y = pd.Series(shuffle_csv["class"])
    y = to_one_hot(y, n_class)
    print(y.shape)
#     print(type(x))
#     print(type(y))
#     print(type(ref))

    return x, ref, y



# In[5]:


def to_one_hot(y, n_class):
    return np.eye(n_class)[y.astype(int)]


# In[10]:


def split_ref_dataset(x_test, y_test, ref_test, dev_ratio):

    test_size = len(x_test)
    print(test_size)
    dev_size = (int)(test_size * dev_ratio)
    print(dev_size)

    x_dev = x_test[:dev_size]
    x_test = x_test[dev_size:]
    y_dev = y_test[:dev_size]
    y_test = y_test[dev_size:]
    ref_dev = ref_test[:dev_size]
    ref_test = ref_test[dev_size:]
    return x_test, x_dev, y_test, y_dev, ref_test, ref_dev, dev_size, test_size - dev_size


# In[7]:


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
    


# In[8]:


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
#         X = torch.matmul(probs.transpose(1,2), V_lookup)
#         X = probs * V_lookup
        X = (probs * V_lookup).sum(dim=1)

        # Right now we can just ignore the fact that we're doing a linear-transform.
        # In the future we'll add nonlinearities

        # Return the logits for the classifier
        return self.W(X)
    


# In[11]:


# load data
x_train, refs_train, y_train = load_ref_data("./data/ref-train-BRAF.csv", sample_ratio=1)
x_test, refs_test, y_test = load_ref_data("./data/ref-test-BRAF.csv", sample_ratio=1)

# split dataset to test and dev
x_train, x_softval, y_train, y_softval, refs_train, refs_softval, softval_size, train_size =     split_ref_dataset(x_train, y_train, refs_train, 0.01)
    
print("Soft Validation size: ", softval_size)
print("Training size: ", train_size)


# In[12]:


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
                                           # shuffle=True, collate_fn=lambda batch: torch.LongTensor(batch))

hardval_loader = torch.utils.data.DataLoader(dataset=hardval_dataset,
                                             batch_size=BATCH_SIZE,
                                             shuffle=True)
                                             # shuffle=True, collate_fn=lambda batch: torch.LongTensor(batch))

softval_loader = torch.utils.data.DataLoader(dataset=softval_dataset,
                                             batch_size=BATCH_SIZE,
                                             shuffle=True)
                                             # shuffle=True, collate_fn=lambda batch: torch.LongTensor(batch))


# In[13]:


# model = AttentionLR(MAX_SEQ_LENGTH, MAX_LABEL)
model = SequenceAttentionClassifier(genome_length=GENOME_LENGTH,
                                    read_length=READ_LENGTH,
                                    vocab_size=VOCAB_SIZE,
                                    query_size=QUERY_SIZE,
                                    embedding_size=EMBEDDING_SIZE,
                                    num_classes=NUM_CLASSES)

criterion = nn.CrossEntropyLoss()
# criterion = nn.L1Loss
optimizer = torch.optim.SGD(model.parameters(), lr=LEARNING_RATE)


# In[15]:


num_epochs = 1

# Training process

b = 0 # count batch
for epoch in range(num_epochs):
    for x_batch, y_batch in train_loader:

        # print(y_batch)
        # print(torch.max(y_batch, 1))
        # print(torch.max(y_batch.type(torch.LongTensor), 1)[1])
        
        optimizer.zero_grad()
        
        outputs = model(x_batch) 

        # loss = criterion(torch.LongTensor(outputs), torch.LongTensor(y_batch))  
        loss = criterion(outputs, torch.max(y_batch.type(torch.LongTensor), 1)[1])
        # loss = criterion(outputs, y_batch)
        loss.backward()
        optimizer.step()
        
        if (b + 1) % 10 == 0:
            print("Epoch {}, Batch {}, loss :{}".format(epoch + 1, b + 1, loss.data[0]))
        b = b + 1
            


# In[ ]:


# # Hard Validation process

# correct = 0
# total = 0
# for i, (x_batch, y_batch) in hardval_loader:
#     outputs = model(x_batch)
#     _, predicted = torch.max(outputs.data, 1)
#     total += labels.size(0)
#     correct += (predicted == labels).sum()

# print('Hard validation Accuracy: {}%'.format(100 * correct / total))    


# # Soft Validation process

# correct = 0
# total = 0
# for i, (x_batch, y_batch) in softval_loader:
#     outputs = model(x_batch)
#     _, predicted = torch.max(outputs.data, 1)
#     total += labels.size(0)
#     correct += (predicted == labels).sum()

# print('Soft validation Accuracy: {}%'.format(100 * correct / total))    

