import sys
import random
from Bio import SeqIO
import numpy as np
import pandas as pd
import sklearn as sk
import tensorflow as tf
import keras
from keras.models import Sequential, Model, load_model
from keras.optimizers import Adam
from keras.callbacks import EarlyStopping, History, ModelCheckpoint
from scipy import stats
from sklearn.metrics import mean_squared_error
import shap

np.random.seed(23)

#Load csv file created from preprocessing and any extra training samples
df = pd.read_csv('preproc_df.csv')
extradf=pd.read_csv('extraPeaks_df.csv')
df = pd.concat([df,extradf])

df['seq']=df['seq'].str.upper()

def train_valid_test_split(dat, valid_frac=0.2, test_frac=0.2,tol=0.1):
    """
    Function to create a train validatin test split between chromosomes.
    Takes all training regions and splits to give proportions based on
    user selected.
    """
    #get the coverage filtered dataset to group by chrom and get num entries
    counts_chroms = dat.groupby('chr')['score'].count().reset_index()
    counts_chroms['prop'] = counts_chroms['score']/counts_chroms['score'].sum()
    def sample_set(counts_chroms,frac = 0.2,tol=0.15):
        #make a copy so can reset
        counts_chroms_ = counts_chroms.copy()
        reset_count=0
        act_frac=0
        sample_chr=[]
        while abs(act_frac-frac)>tol:
            #if gone to high, start again
            if act_frac>frac:
                sample_chr=[]
                act_frac=0
                counts_chroms_ = counts_chroms.copy()
                reset_count+=1
                assert reset_count<100, f"Can't find a combination for {frac} proportion!"
            else:
                #sample rand chrom
                chr_rand = random.randrange(counts_chroms_.shape[0])
                sample_chr.append(counts_chroms_.loc[chr_rand]['chr'])
                act_frac+=counts_chroms_.loc[chr_rand]['prop']
                counts_chroms_ = counts_chroms_.drop(index=[chr_rand]).reset_index(drop=True)
        return(counts_chroms_,sample_chr,act_frac)
    #get test set
    counts_chroms,test_chr,test_act = sample_set(counts_chroms,frac = test_frac,tol=tol)
    #get valid set
    counts_chroms,valid_chr,valid_act = sample_set(counts_chroms,frac = valid_frac,tol=tol)
    #train - remainder
    train_chr = counts_chroms['chr'].tolist()
    train_act = counts_chroms['prop'].sum()
    
    return({"train":train_chr,"valid":valid_chr,"test":test_chr},
           {"train":train_act,"valid":valid_act,"test":test_act}
           )

random.seed(23)

#Train test split
split_dic=train_valid_test_split(df)[0]

trainset=df[df['chr'].isin(split_dic['train'])]
testset=df[df['chr'].isin(split_dic['test'])]
validset=df[df['chr'].isin(split_dic['valid'])]

X_train,Y_train=trainset['seq'],trainset.iloc[:,2:3]
X_test,Y_test=testset['seq'],testset.iloc[:,2:3]
X_valid,Y_valid=validset['seq'],validset.iloc[:,2:3]

#Normalize scores
scaler = MinMaxScaler()
Y_train=scaler.fit_transform(Y_train)
Y_test=scaler.transform(Y_test)
Y_valid=scaler.transform(Y_valid)

#One-hot encoding of nucleotides
nucleotide_dict = {'A': [1, 0, 0, 0],
                   'C': [0, 1, 0, 0],
                   'G': [0, 0, 1, 0],
                   'T': [0, 0, 0, 1],
                   'N': [0, 0, 0, 0]}
def one_hot_encode(seq):
    return np.array([nucleotide_dict[nuc] for nuc in seq])
X_train=np.array(X_train.apply(one_hot_encode).tolist())
X_test=np.array(X_test.apply(one_hot_encode).tolist())
X_valid=np.array(X_valid.apply(one_hot_encode).tolist())

tf.random.set_seed(24)

#Parameters of model
params= {'batch_size': 32, # number of examples per batch
                      'epochs': 100, # number of epochs
                      'early_stop': 20, # patience to reduce training time; you can increase the patience to see if the model improves after more epochs
                      'lr': 0.001, # learning rate
                      'n_conv_layer': 3, # number of convolutional layers
                      'num_filters1': 128, # number of filters/kernels in the first conv layer
                      'num_filters2': 60, # number of filters/kernels in the second conv layer
                      'num_filters3': 60, # number of filters/kernels in the third conv layer
                      'kernel_size1': 7, # size of the filters in the first conv layer
                      'kernel_size2': 3, # size of the filters in the second conv layer
                      'kernel_size3': 5, # size of the filters in the third conv layer
                      'n_dense_layer': 1, # number of dense/fully connected layers
                      'dense_neurons': 64, # number of neurons in the dense layer
                      'dropout_conv': 'yes', # add dropout after convolutional layers?
                      'dropout_prob': 0.4, # dropout probability
                      'pad':'same'}

#Create model
model = Sequential()
model.add(kl.Conv1D(params['num_filters1'], kernel_size=params['kernel_size1'],padding=params['pad'],activation='relu',name='Conv1D_1',input_shape=(251, 4)))
model.add(BatchNormalization())
model.add(MaxPooling1D(2))

model.add(kl.Conv1D(params['num_filters2'], kernel_size=params['kernel_size2'],padding=params['pad'],activation='relu',name='Conv1_2'))
model.add(BatchNormalization())
model.add(MaxPooling1D(2))
model.add(Dropout(params['dropout_prob']))

model.add(kl.Conv1D(params['num_filters3'], kernel_size=params['kernel_size3'],padding=params['pad'],activation='relu',name='Conv1_3'))
model.add(BatchNormalization())
model.add(MaxPooling1D(2))
model.add(Dropout(params['dropout_prob']))

model.add(Flatten())

model.add(Dense(params['dense_neurons'],
                     name=str('Dense_'+str(1)),activation = 'relu'))

model.add(BatchNormalization())
model.add(Dropout(params['dropout_prob']))
model.add(Dense(1, activation='softplus', name=str('Dense_TRansFac')))


model.compile(Adam(learning_rate=params['lr']),loss='mse',metrics=['mse',tf.keras.losses.huber,tf.keras.losses.poisson])

history = model.fit(X_train, Y_train,validation_data=(X_valid, Y_valid),verbose=1, batch_size=params['batch_size'],epochs=params['epochs'],callbacks=[EarlyStopping(patience=params['early_stop'], monitor="val_loss", restore_best_weights=True), History()])

background = X_test[np.random.choice(X_test.shape[0], 1000, replace=False)]

pred = model.predict(X_test, batch_size=params['batch_size'])

print('PCC: ',stats.pearsonr(Y_test.squeeze(), pred.squeeze())[0])

shap.explainers._deep.deep_tf.op_handlers["AddV2"] = shap.explainers._deep.deep_tf.passthrough
explainer = shap.DeepExplainer((model.layers[0].input,model.layers[-1].output),data=background)

shapvalues= explainer.shap_values(X_test)

shapvalues = np.reshape(shapvalues, (shapvalues.shape[0],251,4))
final_contr_scores = shapvalues*X_test

shaplengthlast=np.swapaxes(shapvalues,1,2)
Xlengthlast=np.swapaxes(X_test,1,2)

np.save('shap.npy',shaplengthlast)
np.save('seq.npy', Xlengthlast)
