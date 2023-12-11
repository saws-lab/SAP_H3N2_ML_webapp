# -*- coding: utf-8 -*-
"""
Streamlit-based web application for seasonal antigenic prediction (SAP) of IAV H3N2
using proposed AdaBoost model

Shah et. al., "Seasonal antigenic prediction of influenza A H3N2 using machine
learning", Nature communications.

@author: Awais
"""

#%%
#########
# Imports
#########
import streamlit as st
import pandas as pd
import numpy as np
import joblib

# self defined functions for web-app
import app_utilities


#%%
###########
# variables
###########
mut_mat = 'GIAG010101'   # amino acid mutation matrix identifier
meta_features = [
                 'virus',   # virus avidity
                 'serum',   # antiserum potency
                 'virusPassCat',   # virus passage category
                 'serumPassCat'   # serum passage category
                 ]   # metadata features



#%%
#######################
# Title and information
#######################
st.title('Seasonal antigenic prediction of influenza A virus (IAV) H3N2')
st.write("""
         Predict the NHT-based antigenic difference between virus-antiserum
         pairs using their HA1 sequences and (optional) metadata information.
         ***
         """)


#%%
###########################################
# Select and load trained model and encoder
###########################################

test_seasons = [str(year)+s for year in range (2018, 2021) for s in ['NH', 'SH']] + ['2021NH']

select_season = st.sidebar.selectbox('Select the influenza season of your test virus isolates',
                                     test_seasons,
                                     index=len(test_seasons)-1
                                     )

# as per selected season, find the training end season
# i.e., previous season than the selected season
train_end_season_ind = test_seasons.index(select_season) - 1
train_end_season = test_seasons[train_end_season_ind]

# load the model and encoder as per the train end season
model_fn   = f"trained_model/SAP_H3N2_AdaBoost_model_trained_2003NH_{train_end_season}.joblib"
encoder_fn = f"trained_model/SAP_H3N2_OneHotEncoder_trained_2003NH_{train_end_season}.joblib"

model   = joblib.load(model_fn)
encoder = joblib.load(encoder_fn)


st.sidebar.write("""
                 ***
                 """)



#%%
#######################################
# Sidebar - Single input or file upload
#######################################
st.sidebar.write('''
                 For a single virus-antiserum pair, select 'Input' option,
                 and for multiple pairs, select 'Upload' option.
                 '''
                )

# selectbox with a default option of Input
option = st.sidebar.selectbox('Input data or File upload?',
                              ('Input', 'Upload'),
                              index=0
                             )


# If upload option is selected
if option == 'Upload':
    # get input file and preprocess data
    data = app_utilities.process_file_data(mut_mat)
    
    # proceed only with valid input data
    if data is not None:
        #########################
        # prepare input for model
        #########################
        
        # input features (encoded genetic difference)
        X_test = pd.DataFrame(data.seq_diff.to_list(),
                              index=data.index)
        
        
        #Input features (metadata features)
        X_test_meta  = data[meta_features].fillna('None').astype('str')
        # one hot encoding
        X_test_meta  = encoder.transform(X_test_meta).toarray()
        
        # Combined input
        #    - Encoded genetic difference of HA1 sequences, and
        #    - One-hot encoded metadata information
        # input
        X_test = np.hstack((X_test.values, X_test_meta))
        
        
        ###################
        # get predicted NHT
        ###################
        pred_nht = model.predict(X_test)
        
        data['predicted NHT'] = pred_nht
        
        
        ########################
        # display predicted NHTs
        ########################
        st.subheader('Predicted NHT:')
        st.write(data[['virusName', 'virusPassage', 'serumName', 'serumPassage', 'predicted NHT']])
        
        
        #################
        # Download output
        #################
        out_cols = ['virusName', 'virusPassage', 'serumName', 'serumPassage', 'predicted NHT', 'virusSeq', 'serumSeq']
        output = data[out_cols].to_csv(index=False).encode('utf-8')
        st.download_button("Download",
                           output,
                           "output.csv",
                           "text/csv"
                           )

else:
    # get input and preprocess it
    pair = app_utilities.process_input_data(mut_mat)
    
    # proceed only with valid input data
    if pair is not None:
        #########################
        # prepare input for model
        #########################
        
        # input features (encoded genetic difference)
        # convert to numpy
        X_test = np.array(pair['seq_diff'][0]).reshape(1, -1)
        
        #Input features (metadata features)
        X_test_meta  = pair[meta_features].fillna('None').astype('str')
        # one hot encoding
        X_test_meta  = encoder.transform(X_test_meta).toarray()
        
        # Combined input
        #    - Encoded genetic difference of HA1 sequences, and
        #    - One-hot encoded metadata information
        X_test = np.concatenate((X_test, X_test_meta), axis=1)
        
        
        ###################
        # get predicted NHT
        ###################
        pred_nht = model.predict(X_test)
        
        
        #######################
        # display predicted NHT
        #######################
        st.subheader(f'Predicted NHT: {pred_nht[0]:.2f}')
    



