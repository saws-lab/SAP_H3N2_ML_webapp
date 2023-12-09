# -*- coding: utf-8 -*-
"""
Created on Tue Oct 10 00:27:27 2023

@author: Awais
"""

#%%
'''
Imports
'''
import pandas as pd
import numpy as np
import re
import aaindex
import streamlit as st

# aa index initialization, used for encoding
aaindex.init("./AAindex",index='2')


#%%
def passage_category(passage):
    '''
    Compute passage category based on passage information
    
    Parameters
    ----------
    passage (string): passage of isolate.

    Returns
    -------
    (string): passage category

    '''
    
    '''CELL passage'''
    # contains --> MDCK, M#, MX, IAT, S#, SX, MK, P#, P-#, C#, C #, CX, CELL, X#, RII, MEK 
    # equals   --> C, X
    cell_strings = ['MDCK',
                    r'M\d',
                    'MX',
                    'IAT',
                    r'S\d',
                    'SX',
                    'MK',
                    r'P\d',
                    r'C\d',
                    r'C \d',
                    'CX',
                    'CELL',
                    r'X\d',
                    'RII',
                    'MEK'
                    ]
    
    for check in cell_strings:
        if re.search(check, passage):
            return 'CELL'
    
    if passage == 'C':
        return 'CELL'
    if passage == 'X':
        return 'CELL'
    
    
    '''EGG passage'''
    # To categorize egg passage the search 'E' will not work due to following words
    # passage, material, specimen, plaque, direct, cells, cell, isolate, sample,
    # undefined, vero, uncultered, initial
    #
    # contains --> EX, E#, EGG, SPF 
    # equals   --> E
    
    if re.search(r'E\d', passage):
        return 'EGG'
    if 'EGG' in passage:
        return 'EGG'
    if 'SPF' in passage:
        return 'EGG'
    if passage == 'E':
        return 'EGG'



#%%
def preprocess_input(virus_seq, serum_seq,
                     virus_name=None, serum_name=None,
                     virus_passage=None, serum_passage=None):
    
    '''
    Clean and verify input data
    
    For HA1 sequences
        - remove whitespaces and newlines
        - convert to upper case
        - verify valid amino acids
    
    For name and passage information
        - remove leading and trailing whitespaces
        - convert to upper case
        - get passage category based on passage information
    
    Parameters
    ----------
    virus_seq (string): HA1 sequence of virus isolate.
    serum_seq (string): HA1 sequence of antiserum isolate.
    virus_name (string): Name of virus isolate, default=None.
    serum_name (string): Name of antiserum, default=None.
    virus_passage (string): Passage of virus isolate, default=None.
    serum_passage (string): Passage of antiserum, default=None.

    Returns
    -------
    (dataframe): processed parameters
    '''
    
    # 20 valid amino acids and unknown amino acid 'X'
    aa = ['A','R','N','D','C','Q','E','G','H','I',
          'L','K','M','F','P','S','T','W','Y','V',
          'X']
    
    # to collect processed input data
    pair = {}
    
    # HA1 sequences
    # remove space and newlines
    # convert to upper case letters
    pair['virusSeq'] = virus_seq.replace(" ", "").replace("\n", "").upper()
    pair['serumSeq'] = serum_seq.replace(" ", "").replace("\n", "").upper()
    
    
    # verify virus sequence length should be 329
    if len(pair['virusSeq']) != 329:
        st.write("Invalid length of HA1 sequence of virus isolate.")
        return None
    
    
    # verify antiserum sequence length should be 329
    if len(pair['serumSeq']) != 329:
        st.write("Invalid length of HA1 sequence of antiserum.")
        return None

    
    # verify valid amino acids plus unknown amino acid 'X'
    # for virus sequence
    if not all(x in aa for x in pair['virusSeq']):
        st.write("Invalid amino acids in HA1 sequence of virus isolate.")
        return None
    
    
    # for antiserum sequence
    if not all(x in aa for x in pair['serumSeq']):
        st.write("Invalid amino acids in HA1 sequence of antiserum.")
        return None
    
    
    # name and passage of isolates
    # remove leading and trailing whitespaces
    # convert to upper case
    if virus_name is not None:
        pair['virusName'] = virus_name.strip().upper()
    
    
    if virus_passage is not None:
        pair['virusPassage'] = virus_passage.strip().upper()
        
        # get passage category of virus isolate
        pair['virusPassCat'] = passage_category(pair['virusPassage'])
    
    
    if (virus_name is not None) and (virus_passage is not None):
        # combine name and passage with string '='
        pair['virus'] = pair['virusName'] + '=' + pair['virusPassage']
            
            
    
    # for antiserum
    if serum_name is not None:
        pair['serumName'] = serum_name.strip().upper()
    
    
    if serum_passage is not None:
        pair['serumPassage'] = serum_passage.strip().upper()
        
        # get passage category of antiserum
        pair['serumPassCat'] = passage_category(pair['serumPassage'])
    
    
    if (serum_name is not None) and (serum_passage is not None):
        # combine name and passage with string '='
        pair['serum'] = pair['serumName'] + '=' + pair['serumPassage']
            
    
    # convert to dataframe
    pair = pd.DataFrame.from_dict([pair])
    
    return pair

   
 
#%%
def seq_encode(virusSeq, serumSeq, mut_mat):
    '''
    Compute pairwise genetic difference at each site of virusSeq and serumSeq as per aaIndID
    return a list of same length as of sequence, where each element of list is
    a difference at each site of virus and serum sequences
    
    Note: Here, we assume that valid amino acid characters are already verified
          in virusSeq and serumSeq
    
    Parameters
    ----------
    virusSeq (string): amino acid sequence of virus isolate
    serumSeq (string): amino acid sequence of antiserum
    mut_mat (string): amino acid mutation matrix identifier
    
    Returns
    -------
    diff (list): difference at each site of virusSeq and serumSeq
    '''
    
    # list of differences at each site of virusSeq and serumSeq
    diff = []
    
    # get the similarity matrix as per AA Index ID
    Sim_Mat = aaindex.get(mut_mat)
    
    # loop through each site of a sequence
    for i in range(len(virusSeq)):
        # unknown character 'X' encoded as 0
        if virusSeq[i]=='X' or serumSeq[i]=='X':
            tmp = 0
        # for remaining valid amino acid characters
        else:
            # Encoding as per mut_mat similarity matrix
            m = Sim_Mat.get(virusSeq[i], virusSeq[i])
            n = Sim_Mat.get(serumSeq[i], serumSeq[i])
            k = Sim_Mat.get(serumSeq[i], virusSeq[i])
    
            tmp = m + n - 2*k      # difference/dissimilarity formula

        # save distance at site i
        diff.append(tmp)
        
    return diff



#%%
def process_input_data(mut_mat):
    '''
    Get single virus-antiserum pair and prepare it for model input
    
    Parameters
    ----------
    mut_mat (string): amino acid mutation matrix identifier
    
    Returns
    -------
    pair (dataframe): processed input for model
       
    '''
    
    '''
    Header
    '''
    st.header('Input HA1 sequences and (optional) metadata')


    '''
    Example inputs
    '''
    tmp_virus_seq = 'QKIPGNDNSTATLCLGHHAVPNGTIVKTITNDRIEVTNATELVQNSSIGEICDSPHQILDGGNCTLIDALLGDPQCDGFQNKEWDLFVERSRANSNCYPYDVPDYASLRSLVASSGTLEFKNESFNWTGVKQNGTSSACIRGSSSSFFSRLNWLTHLNNIYPAQNVTMPNKEQFDKLYIWGVHHPDTDKNQISLFAQSSGIITVSTKRSQQAVIPNIGSRPRIRDIPSRISIYWTIVKPGDILLINSTGNLIAPRGYFKIRSGKSSIMRSDAPIGKCKSECITPNGSIPNDKPFQNVNRITYGACPRYVKQSTLKLATGMRNVPEKQTR'
    tmp_serum_seq = 'QKIPGNDNSTATLCLGHHAVPNGTIVKTITNDRIEVTNATELVQNSSIGEICDSPHQILDGENCTLIDALLGDPQCDGFQNKKWDLFVERNKAYSNCYPYDVPDYASLRSLVASSGTLEFNNESFNWAGVTQNGTSSSCIRGSKSSFFSRLNWLTHLNSKYPALNVTMPNNEQFDKLYIWGVHHPGTDKNQISLYAQSSGRITVSTKRSQQAVIPNIGSRPRIRDIPSRISIYWTIVKPGDILLIXSTGNLIAPRGYFKIRSGKSSIMRSDAPIGKCKSECITPNGSIPNDKPFQNVNRITYGACPRYVKQSTLKLATGMRNVPERQTR'
    tmp_virus_name = 'A/MICHIGAN/173/2020'
    tmp_serum_name = 'A/KANSAS/14/2017'
    tmp_virus_passage = 'S1'
    tmp_serum_passage = 'E5'


    '''
    Take inputs from user
        - virus HA1 sequence
        - antiserum HA1 sequence
        - name of virus isolate (optional)
        - name of antiserum (optional)
        - passage of virus isolate (optional)
        - passage of antiserum (optional)
    '''
    st.text_area('HA1 sequence of virus isolate',
                 tmp_virus_seq,
                 height=140,
                 max_chars=329,
                 key='virus_seq')
    st.text_area('HA1 sequence of antiserum/vaccine (reference virus isolate)',
                 tmp_serum_seq,
                 height=140,
                 max_chars=329,
                 key='serum_seq')
    st.text_input('Name of virus isolate (if unknown, leave it blank)',
                  tmp_virus_name,
                  key='virus_name')
    st.text_input('Name of antiserum/vaccine (if unknown, leave it blank)',
                  tmp_serum_name,
                  key='serum_name')
    st.text_input('Passage of virus isolate (if unknown, leave it blank)',
                  tmp_virus_passage,
                  key='virus_passage')
    st.text_input('Passage of antiserum/vaccine (if unknown, leave it blank)',
                  tmp_serum_passage,
                  key='serum_passage')


    # prediction button
    if st.button('Predict',
                 help='Predict NHT-based antigenic difference of virus-antiserum pair',
                ):
        # preprocess input
        pair = preprocess_input(st.session_state.virus_seq,
                                st.session_state.serum_seq,
                                virus_name=st.session_state.virus_name,
                                serum_name=st.session_state.serum_name,
                                virus_passage=st.session_state.virus_passage,
                                serum_passage=st.session_state.serum_passage
                                )
        
        # if no input format error occured
        # procees with valid input
        if pair is not None:
            # encode input HA1 sequences using mutation matrix
            seq_diff = seq_encode(pair.virusSeq.values[0],
                                  pair.serumSeq.values[0],
                                  mut_mat
                                  )
            
            pair['seq_diff'] = [seq_diff]
        
            return pair
        
        # in case of input error
        else:
            return None
    
    # When predict button is not used
    else:
        return None



#%%
def process_file_data(mut_mat):
    '''
    Get data of multiple virus-antiserum pairs and prepare it for model input
    
    Parameters
    ----------
    mut_mat (string): amino acid mutation matrix identifier
    
    Returns
    -------
    data_clean (dataframe): processed input data for model
    
    '''
    
    '''
    File instructions
    '''
    st.header('Upload HA1 sequences and (optional) metadata information of virus-antiserum pairs')
    
    # Instructions for file formatting
    st.write('''
             The file should be in **CSV** format with the following headers and
             each row should represent the input data of a single virus-antiserum pair:
             - **virusSeq**: HA1 sequence of virus isolate
             - **serumSeq**: HA1 sequence of antiserum antiserum/vaccine (reference virus isolate)
             - **virusName**: (optional) Name of virus isolate (if unknown, leave it blank)
             - **serumName**: (optional) Name of antiserum/vaccine (if unknown, leave it blank)
             - **virusPassage**: (optional) Passage of virus isolate (if unknown, leave it blank)
             - **serumPassage**: (optional) Passage of antiserum/vaccine (if unknown, leave it blank)
             ''')
    
    # get file from user
    uploaded_file = st.file_uploader("**Upload file**", type=['csv'])
    
    ####################
    # with uploaded file
    ####################
    if uploaded_file is not None:
        # read file into dataframe
        data = pd.read_csv(uploaded_file)
        
        # message for successfully uploaded file
        st.success('File successfully uploaded!', icon="✅")
        
        # proceed if mandatory headers for sequences are found
        if ('virusSeq' in data.columns) and ('serumSeq' in data.columns):
            #################
            # prepare headers
            #################
            # required headers
            req_headers = ['virusSeq', 'serumSeq',
                           'virusName', 'serumName',
                           'virusPassage', 'serumPassage']
            
            # provided headers in uploaded file
            provided_headers = data.columns.intersection(req_headers)
            # get data of provided headers
            data = data[provided_headers]
            
            # add remaining headers
            rem_headers = list(set(req_headers) - set(provided_headers))
            data.loc[:, rem_headers] = None
            
            
            ###############################
            # for each virus-antiserum pair
            #    - preprocess input
            #    - encode HA1 sequences
            ###############################
            # to collect processed data
            data_clean = pd.DataFrame()
            
            for ind, row in data.iterrows():
                # preprocess input
                pair = preprocess_input(row.virusSeq,
                                        row.serumSeq,
                                        virus_name=row.virusName,
                                        serum_name=row.serumName,
                                        virus_passage=row.virusPassage,
                                        serum_passage=row.serumPassage
                                        )
                
                # proceed only with valid input data
                if pair is not None:
                    # encode input HA1 sequences using mutation matrix
                    seq_diff = seq_encode(pair.virusSeq.values[0],
                                          pair.serumSeq.values[0],
                                          mut_mat
                                          )
                    
                    pair['seq_diff'] = [seq_diff]
                    
                    
                    # combine data
                    data_clean = pd.concat((data_clean, pair), ignore_index=True)
                
                # in case of error
                else:
                    st.error(f'At pair number: {ind+1}')
                    return None
            
            
        ####################################
        # if headers for sequences not found
        ####################################
        else:
            # report error
            st.error('Wrong header(s) for input HA1 sequence(s). Check formatting instructions above.')
            return None
        
    ##################
    # no uploaded file
    ##################
    else:
        return None

    
    return data_clean