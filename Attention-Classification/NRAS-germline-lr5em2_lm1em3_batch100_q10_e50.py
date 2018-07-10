
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


isRetrain = False


# Hyperparameter

QUERY_SIZE = 10
EMBEDDING_SIZE = 50

# HIDDEN_SIZE = 512
# ATTENTION_SIZE = 64

parameter_name = 'lr5em2_lm1em4_batch100_q1_e50'

LEARNING_RATE = 5*1e-2
BATCH_SIZE = 100
LAMBDA = 1e-3

num_epochs = 100

VOCAB_SIZE = 5
NUM_CLASSES = 2

# Data-specific

READ_LENGTH = 100

# # KRAS
# GENOME_START = 25204789
# GENOME_END = 25250936

# # NOTCH1
# GENOME_START = 136494433 
# GENOME_END = 136545786

# NRAS
file_name = 'NRAS-germline'
GENOME_START = 114704464 
GENOME_END = 114716894

GENOME_LENGTH = GENOME_END - GENOME_START + 1  
CONTEXT_SIZE = GENOME_LENGTH

ref_names = ["class", "ref", "sequence"]
VOCAB = ['N','A','T','C','G']


# In[3]:


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


# In[4]:


def create_synthetic_data(file_name, cancer_genes, benign_genes, num_patients=10, num_reads_per_patients=3, read_length=READ_LENGTH, genome_length=GENOME_LENGTH, vocab=VOCAB, isSomatic=True, print_seq=False):
    
    seq_list = np.random.choice(vocab, [num_patients, genome_length], replace=True)
    backup_seq_list = seq_list
    
    for loc, mutation in cancer_genes.items():
        seq_list[np.random.choice(num_patients, int(num_patients*mutation[1]), replace=False), loc] = mutation[0]
    
    genomes = []
    for r in range(seq_list.shape[0]):
        seq = ''.join(seq_list[r,:])
        if print_seq:
            print(seq)
        genomes.append(seq)

    locs = np.random.choice(genome_length-read_length, num_patients*num_reads_per_patients)

    file = open('./tumor-genome-'+file_name+'.txt','w')
    count = 0
    reads = []
    for genome in genomes:
        for t in range(num_reads_per_patients):
            index = count*num_reads_per_patients+t
            reads.append(genome[locs[index]:locs[index]+read_length])
        file.write("%s\n" % genome)
        count = count + 1
    file.close() 

    file = open('./tumor-syn-'+file_name+'.csv','w')
    for r in range(num_patients*num_reads_per_patients):
        file.write("1, %d, %s\n" % (locs[r], reads[r]))
    file.close()
    
    tumor_locs = locs
    tumor_reads = reads

    if isSomatic:
        seq_list = backup_seq_list
    else:
        seq_list = np.random.choice(vocab, [num_patients, genome_length], replace=True)
    
    for loc, mutation in benign_genes.items():
        seq_list[np.random.choice(num_patients, int(num_patients*mutation[1]), replace=False), loc] = mutation[0]
    
    genomes = []
    for r in range(seq_list.shape[0]):
        seq = ''.join(seq_list[r,:])
        if print_seq:
            print(seq)
        genomes.append(seq)

    locs = np.random.choice(genome_length-read_length, num_patients*num_reads_per_patients)
    
    file = open('./normal-genome-'+file_name+'.txt','w')
    count = 0
    reads = []
    for genome in genomes:
        for t in range(num_reads_per_patients):
            index = count*num_reads_per_patients+t
            reads.append(genome[locs[index]:locs[index]+read_length])
        file.write("%s\n" % genome)
        count = count + 1
    file.close() 

    file = open('./normal-syn-'+file_name+'.csv','w')
    for r in range(num_patients*num_reads_per_patients):
        file.write("0, %d, %s\n" % (locs[r], reads[r]))
    file.close() 
    
    normal_locs = locs
    normal_reads = reads
    
    file = open('./syn-'+file_name+'.csv','w')
    for r in range(num_patients*num_reads_per_patients):
        file.write("1,%d,%s\n" % (tumor_locs[r], tumor_reads[r]))
        file.write("0,%d,%s\n" % (normal_locs[r], normal_reads[r]))
    file.close() 
    
    return './syn-'+file_name+'.csv'
    
    
    


# In[5]:


def to_one_hot(y, n_class):
    
    return np.eye(n_class)[y.astype(int)]


# In[6]:


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
        
        # self.read_length = read_length
        self.labels = labels
        self.reads = reads
        self.ref_locs = ref_locs
        self.genome_start = genome_start
        self.genome_end = genome_end

    def __len__(self):
        return len(self.reads)

    def __getitem__(self, index):
        
        vals = list(self.reads[index])
        read_length = len(vals)
        locs = list(np.arange(self.ref_locs[index]-self.genome_start,self.ref_locs[index]+read_length-self.genome_start))

#         print(len(vals))
#         print(len(locs))
        
        vals2idx = {'N': 0, 'A': 1, 'C': 2, 'T': 3, 'G': 4}
#         read = torch.LongTensor(np.array([vals2idx[val]+loc*len(vals2idx) for val, loc in zip(vals, locs)], dtype=int), requires_grad=False)

        read = torch.autograd.Variable(torch.LongTensor(np.array([vals2idx[val]+loc*len(vals2idx) for val, loc in zip(vals, locs)], dtype=int)), requires_grad=False)
        
        X = read
        Y = self.labels[index,:]

#         torch.LongTensor(self.labels[index,:])

        return X, Y
    


# In[8]:


class SequenceAttentionClassifier(nn.Module):
    
    def __init__(self, genome_length, vocab_size=5, query_size=10, embedding_size=128, num_classes=2):
        
        super(SequenceAttentionClassifier, self).__init__()
        
        self.genome_length = genome_length
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
        probs = F.softmax(logits, dim = 0)
        
        # Calculate the covariates for the logistic regression
#         X = torch.matmul(probs.transpose(1,2), V_lookup)
#         X = probs * V_lookup
        X = (probs * V_lookup).sum(dim=1)

        # Right now we can just ignore the fact that we're doing a linear-transform.
        # In the future we'll add nonlinearities

        # Return the logits for the classifier
        return self.W(X), K_lookup, V_lookup
    


# In[9]:


def GetTensorBatch(reads, ref_locs, labels, batch_size=100, genome_start=0, genome_end=0):
      
    batches = {}
    set_size = len(ref_locs)
    for batch in range(set_size // batch_size):
        
        x_batch = []
        y_batch = []
        
        for index in range(batch_size):
        
            vals = list(reads[index])
            read_length = len(vals)
            locs = list(np.arange(ref_locs[index]-genome_start,ref_locs[index]+read_length-genome_start))
        
            vals2idx = {'N': 0, 'A': 1, 'C': 2, 'T': 3, 'G': 4}
            read = torch.autograd.Variable(torch.LongTensor(np.array([vals2idx[val]+loc*len(vals2idx) for val, loc in zip(vals, locs)], dtype=int)), requires_grad=False)
        
            X = read
            Y = labels[index,:]
            x_batch.append(X)
            y_batch.append(Y)
        
        batches[batch] = [x_batch, y_batch]

    return batches
    


# In[14]:


# load data

x_train, refs_train, y_train = load_ref_data("../data/ref-germline-NRAS-train39000.csv", sample_ratio=1)
x_hardtest, refs_hardtest, y_hardtest = load_ref_data("../data/ref-germline-NRAS-test1000.csv", sample_ratio=1)

# split dataset to train and validation
x_train, x_test, y_train, y_test, refs_train, refs_test, test_size, train_size = split_ref_dataset(x_train, y_train, refs_train, 0.2)
    
# split dataset to test and dev
x_softtest, x_val, y_softtest, y_val, refs_softtest, refs_val, val_size, softtest_size = split_ref_dataset(x_test, y_test, refs_test, 0.5)

print("Training size: ", train_size)
print("Soft Test size: ", softtest_size)
print("Hard Test size: ", len(y_hardtest))
print("Validation size: ", val_size)

# In[15]:


# Generators
# train_dataset = TensorizedReadDataset(reads=x_train, 
#                                       ref_locs=refs_train, 
#                                       labels=y_train, 
#                                       read_length=READ_LENGTH, 
#                                       genome_start=GENOME_START, 
#                                       genome_end=GENOME_END)

# hardval_dataset = TensorizedReadDataset(reads=x_test, 
#                                         ref_locs=refs_test, 
#                                         labels=y_test, 
#                                         read_length=READ_LENGTH, 
#                                         genome_start=GENOME_START, 
#                                         genome_end=GENOME_END)

# softval_dataset = TensorizedReadDataset(reads=x_softval, 
#                                         ref_locs=refs_softval, 
#                                         labels=y_softval, 
#                                         read_length=READ_LENGTH, 
#                                         genome_start=GENOME_START, 
#                                         genome_end=GENOME_END)

# # Input pipeline
# train_loader = torch.utils.data.DataLoader(dataset=train_dataset,
#                                            batch_size=BATCH_SIZE,
#                                            shuffle=True)

# hardtest_loader = torch.utils.data.DataLoader(dataset=hardval_dataset,
#                                              batch_size=BATCH_SIZE,
#                                              shuffle=True)

# softval_loader = torch.utils.data.DataLoader(dataset=softval_dataset,
#                                              batch_size=BATCH_SIZE,
#                                              shuffle=True)

train_loader = GetTensorBatch(reads=x_train, 
                                      ref_locs=refs_train, 
                                      labels=y_train, 
                                      batch_size=BATCH_SIZE,
                                      genome_start=GENOME_START, 
                                      genome_end=GENOME_END)

val_loader = GetTensorBatch(reads=x_val, 
                                      ref_locs=refs_val, 
                                      labels=y_val, 
                                      batch_size=BATCH_SIZE,
                                      genome_start=GENOME_START,
                                      genome_end=GENOME_END)

hardtest_loader = GetTensorBatch(reads=x_hardtest, 
                                        ref_locs=refs_hardtest, 
                                        labels=y_hardtest, 
                                        batch_size=BATCH_SIZE, 
                                        genome_start=GENOME_START, 
                                        genome_end=GENOME_END)

softtest_loader = GetTensorBatch(reads=x_softtest, 
                                        ref_locs=refs_softtest, 
                                        labels=y_softtest, 
                                        batch_size=BATCH_SIZE,
                                        genome_start=GENOME_START, 
                                        genome_end=GENOME_END)



# In[16]:


# isRetrain = True


# In[17]:


model = SequenceAttentionClassifier(genome_length=GENOME_LENGTH,

                                    vocab_size=VOCAB_SIZE,
                                    query_size=QUERY_SIZE,
                                    embedding_size=EMBEDDING_SIZE,
                                    num_classes=NUM_CLASSES)

criterion = nn.CrossEntropyLoss()
optimizer = torch.optim.SGD(model.parameters(), lr=LEARNING_RATE)

if isRetrain:
    model.load_state_dict(torch.load('./'+file_name+'_model.pth'))


# In[18]:

f1 = open('cv_acc_each_epoch_'+parameter_name+'.txt', 'a+')
f2 = open('cv_loss_each_batch_'+parameter_name+'.txt', 'a+')
f3 = open('cv_loss_val_each_batch_'+parameter_name+'.txt', 'a+')
f4 = open('cv_acc_each_batch_'+parameter_name+'.txt', 'a+')

lam = LAMBDA

# Training process
for epoch in range(num_epochs):
#     b = 0 # count batch
    for b, [x_batch, y_batch] in train_loader.items():
        
        full_loss = 0
        full_data_loss = 0
        optimizer.zero_grad()
        for c in range(BATCH_SIZE):
            
            x_input = x_batch[c].view(1, x_batch[c].shape[0])
            y_input = torch.Tensor(y_batch[c]).type(torch.float64).view(1, y_batch[c].shape[0])
        
            outputs, K, V = model(x_input)
            data_loss = criterion(outputs, torch.max(y_input, 1)[1])
            
            reg_loss = (torch.norm(K,2,1).mean() + torch.norm(V,2,1).mean())*lam
            loss = data_loss + reg_loss
            full_loss = full_loss + loss
            full_data_loss = full_data_loss + data_loss

            loss.backward()
            optimizer.step()
            f2.write(repr(data_loss)+"\n")
        
#         full_loss = full_loss / 2
        # full_loss.backward()
        # optimizer.step()
        
        # if (b + 1) % 1 == 0:
        print("Epoch {}, Batch {}, loss :{}".format(epoch + 1, b + 1, full_loss))
#         b = b + 1
        
        # loss_each_batch.append(full_loss)
        # f2.write(repr(full_data_loss)+"\n")

        full_loss = 0
        full_data_loss = 0
        
        correct = 0
        total = 0
        for b, [x_batch, y_batch] in val_loader.items():
            full_vloss = 0
            for c in range(BATCH_SIZE):
                
                x_input = x_batch[c].view(1, x_batch[c].shape[0])
                y_input = torch.Tensor(y_batch[c]).type(torch.float64).view(1, y_batch[c].shape[0])
                
                outputs, K, V = model(x_input)
                _, predicted = torch.max(outputs.data, 1)
                
                data_vloss = criterion(outputs, torch.max(y_input, 1)[1])
                full_vloss = full_vloss + data_vloss

                total += len(y_input)
                correct += (predicted == torch.max(y_input.type(torch.LongTensor), 1)[1]).sum()

                f3.write(repr(data_vloss)+"\n")
            # f3.write(repr(full_vloss)+"\n")
                    
        acc = int(correct) / total
        print('Validation Accuracy: {}%'.format(100 * acc))   
        # acc_each_epoch.append(acc)
        f4.write(repr(acc)+"\n")
    
    correct = 0
    total = 0
    for b, [x_batch, y_batch] in val_loader.items():
        for c in range(BATCH_SIZE):
            
            x_input = x_batch[c].view(1, x_batch[c].shape[0])
            y_input = torch.Tensor(y_batch[c]).type(torch.float64).view(1, y_batch[c].shape[0])
            
            outputs, K, V = model(x_input)
            _, predicted = torch.max(outputs.data, 1)
            total += len(y_input)
            correct += (predicted == torch.max(y_input.type(torch.LongTensor), 1)[1]).sum()
    acc = int(correct) / total
    print('Validation Accuracy: {}%'.format(100 * acc))   
    # acc_each_epoch.append(acc)
    f1.write(repr(acc)+"\n")

f1.close() 
f2.close() 
f3.close() 
f4.close() 

# In[202]:

torch.save(model, './full_' + file_name + '_' + parameter_name + '_model.pt')
torch.save(model.state_dict(), './'+file_name + '_' + parameter_name + '_model.pth')



# In[ ]:


# plt.rcParams['figure.figsize'] = [30, 5]

# plt.plot(np.arange(len(loss_each_batch))+1, loss_each_batch)
# plt.xlim(1, len(loss_each_batch))
# plt.title('loss vs. batchs')
# plt.show()

# plt.plot(np.arange(len(acc_each_batch))+1, acc_each_batch)
# plt.xlim(1, len(acc_each_batch))
# plt.title('accuracy vs. batchs')
# plt.show()

# plt.plot(np.arange(len(acc_each_epoch))+1, acc_each_epoch)
# plt.xlim(1, len(acc_each_epoch))
# plt.title('accuracy vs. epochs')
# plt.show()


# In[10]:


# model = torch.load('./full_'+file_name+'_model.pt')
# model.load_state_dict(torch.load('./'+file_name+'_model.pth'))


# In[16]:



correct = 0
total = 0
for b, [x_batch, y_batch] in softtest_loader.items():
    for c in range(BATCH_SIZE):
    
        x_input = x_batch[c].view(1, x_batch[c].shape[0])
        y_input = torch.Tensor(y_batch[c]).type(torch.float64).view(1, y_batch[c].shape[0])
            
        outputs, K, V = model(x_input)
        _, predicted = torch.max(outputs.data, 1)
        total += len(y_input)
        correct += (predicted == torch.max(y_input.type(torch.LongTensor), 1)[1]).sum()
acc = int(correct) / total
print('Soft Test Accuracy: {}%'.format(100 * acc))   


correct = 0
total = 0
for b, [x_batch, y_batch] in hardtest_loader.items():
    for c in range(BATCH_SIZE):
        
        x_input = x_batch[c].view(1, x_batch[c].shape[0])
        y_input = torch.Tensor(y_batch[c]).type(torch.float64).view(1, y_batch[c].shape[0])
            
        outputs, K, V = model(x_input)
        _, predicted = torch.max(outputs.data, 1)
        total += len(y_input)
        correct += (predicted == torch.max(y_input.type(torch.LongTensor), 1)[1]).sum()
acc = int(correct) / total
print('Hard Test Accuracy: {}%'.format(100 * acc))

