import matplotlib.pyplot as plt

rootf = '.' ## Your folder ##

import sys
sys.path.append(rootf + '/diffRBM/source/')
sys.path.append(rootf + '/diffRBM/utilities/')

import argparse
parser = argparse.ArgumentParser() 
parser.add_argument('-pep', type = str, required=True)
parser.add_argument('-train', type = int, default = 0, required=False, help='Whether to train a model with the training data')
parser.add_argument('-alpha', type = int, default = 0, required=False, help='Choose CDR3 on alpha chain instead of beta chain')
parser.add_argument('-LR', type = int, default = 0, required=False, help='Choose the Left-Right format instead of the alignment')
parser.add_argument('-hu', nargs = '+', type = int, required=False,  default = [20], help='Number of Top RBM hidden units')
parser.add_argument('-l12', nargs = '+', type = float, required=False, default = [0.01], help='L^1_2 Top RBM regularization')
parser.add_argument('-verbose', type = int, required=False, default = 1, help='Switch off verbose mode')
args = parser.parse_args()

makeVJ = 1  ## just to decide if to include V and J or not ##
maketrainingB=0
makeLR = 0
alpha_chain = 0
vv = 1
maketraining=0
makeRW=0

if args.alpha:
    alpha_chain = 1
if args.LR:
    makeLR = 1
    
if args.train:
    maketraining=1
    
if args.verbose==0:
    vv=0


listreg = list(args.l12) 
listhu = list(args.hu) 

pep=args.pep

## Functions ##

import numpy as np
import rbm, diffrbm, utilities, sequence_logo, Proteins_utils, RBM_utils
import pandas as pd
import os
import random

curr_int = np.int16

def produce_callback_weights(RBMdiff):
    def callback():
        RBMdiff.RBMpost.weights[RBMdiff.n_h_:] = 0
    return callback

def gene_to_num_str(gene_name, gene_type):
    """Strips excess gene name info to number string.

    Parameters
    ----------
    gene_name : str
        Gene or allele name
    gene_type : char
        Genomic cassette type. (i.e. V, D, or J)
    Returns
    -------
    num_str : str
        Reduced gene or allele name with leading zeros and excess
        characters removed.
        
    Taken from SONIA: https://github.com/statbiophys/SONIA

    """
    # get rid of allele
    gene_name=gene_name.split('*')[0]
    num_str = gene_type.lower().join([g.lstrip('0') for g in gene_name.lower().split(gene_type.lower())[1:]])
    num_str = '-'.join([g.lstrip('0') for g in num_str.split('-')])
    return gene_type.lower() + num_str.replace('/', '')

def add_VJ_info_2num(seqs_2num, V_list, J_list, V_dict, J_dict):
    Vlist_2num = np.array([V_dict[v] for v in V_list]).reshape((seqs_2num.shape[0], 1))
    Jlist_2num = np.array([J_dict[j] for j in J_list]).reshape((seqs_2num.shape[0], 1))
    seqs_2num_VJ = np.append(seqs_2num, Vlist_2num, axis=1)
    seqs_2num_VJ = np.append(seqs_2num_VJ, Jlist_2num, axis=1)
    return seqs_2num_VJ


## Vtypes and Jtypes for Left-Right format ##
vtypes = ['v10-1', 'v10-2', 'v10-3', 'v11-1', 'v11-2', 'v11-3', 'v12-3', 'v12-4', 'v12-5', 'v13', 'v14', 'v15', 'v16', 'v18', 'v19', 'v2', 'v20-1', 'v24-1', 'v25-1', 'v27', 'v28', 'v29-1', 'v3-1', 'v30', 'v4-1', 'v4-2', 'v4-3', 'v5-1', 'v5-4', 'v5-5', 'v5-6', 'v5-8', 'v6-1', 'v6-2', 'v6-3', 'v6-4', 'v6-5', 'v6-6', 'v6-8', 'v6-9', 'v7-2', 'v7-3', 'v7-4', 'v7-6', 'v7-7', 'v7-8', 'v7-9', 'v9']
jtypes = ['j1-1', 'j1-2', 'j1-3', 'j1-4', 'j1-5', 'j1-6', 'j2-1', 'j2-2', 'j2-3', 'j2-4', 'j2-5', 'j2-6', 'j2-7']
ltypes = ['l' + str(n) for n in range(5,31)]


def encodeLR(full_table, max_L, withVJ=False):
    data_seqs = full_table

    ## Left-Right version ##
    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
    aadict = {amino_acids[k]: k for k in range(len(amino_acids))} 
    aadict_inv = {k: amino_acids[k] for k in range(len(amino_acids))}
    max_depth = max_L
    
    if not(withVJ):

        posl=-1
        posv=100
        posj=100

        seq_train = list(np.array(data_seqs)[:,0])
        seq_train_LR = []
        lenn=[]
        ind_prob=[]
        for t in range(len(seq_train)):
            seq = seq_train[t]
            seq_feature_lsts = [['a' + aa + str(i)] for i, aa in enumerate(seq)]
            seq_feature_lsts += [['a' + aa + str(-1-i)] for i, aa in enumerate(seq[::-1])]
            seq_feature_lsts += [['l' + str(len(seq))]]

            length_input = 2*max_L
            data_enc = np.zeros(length_input, dtype=np.int8)
            positions = [int(seq_feature_lsts[i][0][2:]) for i in range(len(seq_feature_lsts[:posl]))]
            #max_depth1 = int(len(seq)/2)
            #max_depth2 = len(seq) - max_depth1
            c=0
            for i in list(np.arange(max_depth)) + list(np.arange(-max_depth, 0)):
                #if (i in positions and i in list(np.arange(max_depth1))) or (i in positions and i in list(np.arange(-max_depth2,0))):
                if i in positions:
                    ii =  positions.index(i)
                    data_enc[i] = aadict[seq_feature_lsts[ii][0][1]]
                    c+=1
                else:
                    data_enc[i] = len(amino_acids)
                    
            if c!=len(seq[0]):
                ind_prob.append(t)
            seq_train_LR.append(list(data_enc))
            lenn.append(seq_feature_lsts[posl])
        ltypes0 = list(np.unique([lenn[i] for i in range(len(lenn))]))
        
        for v in ltypes0:
            if v not in ltypes:
                print('Unrecognized ltype!')
        NL = len(ltypes)

        def convert_numberL(seqs): # convert to numbers the L types
            ldict = {ltypes[k]: 5+k for k in range(len(ltypes))}
            l_num = np.array(list(map(lambda x: [ldict[y] for y in x], seqs[0:])), dtype=curr_int, order="c")
            return l_num

        seq_train_n = np.array(seq_train_LR)
        lsegs_n = convert_numberL(lenn)

        final_set = []
        for f in range(len(seq_train_n)):
            final_set.append(list(np.hstack((seq_train_n[f],lsegs_n[f]))))
        final_set = np.array(final_set)

    if withVJ:

        posl=-3
        posv=-2
        posj=-1

        seq_train_n=[]
        vsegs=[]
        jsegs=[]
        lenn=[]
        ind_prob=[]
        for t in range(len(data_seqs)):
            seq = data_seqs[t]
            v_genes = [gene.split('*')[0] for gene in seq[1:] if 'v' in gene.lower()]
            j_genes = [gene.split('*')[0] for gene in seq[1:] if 'j' in gene.lower()]

            seq_feature_lsts = [['a' + aa + str(i)] for i, aa in enumerate(seq[0])]
            seq_feature_lsts += [['a' + aa + str(-1-i)] for i, aa in enumerate(seq[0][::-1])]
            seq_feature_lsts += [['l' + str(len(seq[0]))]]
            v_genes = [gene.split('*')[0] for gene in seq[1:] if 'v' in gene.lower()]
            j_genes = [gene.split('*')[0] for gene in seq[1:] if 'j' in gene.lower()]
            try:
                seq_feature_lsts += [[gene_to_num_str(gene,'V')] for gene in v_genes]
                seq_feature_lsts += [[gene_to_num_str(gene,'J')] for gene in j_genes]
            except ValueError:
                pass

            length_input = 2*max_L
            data_enc = np.zeros(length_input, dtype=np.int8)
            positions = [int(seq_feature_lsts[i][0][2:]) for i in range(len(seq_feature_lsts[:posl]))]
            #max_depth1 = int(len(seq[0])/2)
            #max_depth2 = len(seq[0]) - max_depth1
            c=0
            for i in list(np.arange(max_depth)) + list(np.arange(-max_depth, 0)):
                #if (i in positions and i in list(np.arange(max_depth1))) or (i in positions and i in list(np.arange(-max_depth2,0))):
                if i in positions:
                    ii =  positions.index(i)
                    data_enc[i] = aadict[seq_feature_lsts[ii][0][1]]
                    c+=1
                else:
                    data_enc[i] = len(amino_acids)
            #if c!=len(seq[0]):
            #    ind_prob.append(t)
                
            seq_train_n.append(list(data_enc))

            vsegs.append(seq_feature_lsts[posv])
            jsegs.append(seq_feature_lsts[posj])
            lenn.append(seq_feature_lsts[posl])

        vtypes0 = list(np.unique([vsegs[i][0] for i in range(len(vsegs))])) # fixed
        jtypes0 = list(np.unique([jsegs[i][0] for i in range(len(jsegs))]))
        ltypes0 = list(np.unique([lenn[i] for i in range(len(lenn))]))
        
        for v in ltypes0:
            if v not in ltypes:
                print('Unrecognized ltype!')
        for v in vtypes0:
            if v not in vtypes:
                print('Unrecognized vtype!')
        for j in jtypes0:
            if j not in jtypes:
                print('Unrecognized jtype!')

        NJ = len(jtypes)
        NV = len(vtypes)
        NL = len(ltypes)

        def convert_numberV(seqs): # convert to numbers the V types, where the 'types' alphabet is sample-dependent
            vdict = {vtypes[k]: k for k in range(len(vtypes))}
            v_num = np.array(list(map(lambda x: [vdict[y] for y in x], seqs[0:])), dtype=curr_int, order="c")
            return v_num

        def convert_numberJ(seqs): # convert to numbers the J types
            jdict = {jtypes[k]: k for k in range(len(jtypes))}
            j_num = np.array(list(map(lambda x: [jdict[y] for y in x], seqs[0:])), dtype=curr_int, order="c")
            return j_num

        def convert_numberL(seqs): # convert to numbers the L types
            ldict = {ltypes[k]: 5+k for k in range(len(ltypes))}
            l_num = np.array(list(map(lambda x: [ldict[y] for y in x], seqs[0:])), dtype=curr_int, order="c")
            return l_num

        seq_train_n = np.array(seq_train_n)
        vsegs_n = convert_numberV(vsegs)
        jsegs_n = convert_numberJ(jsegs)
        lsegs_n = convert_numberL(lenn)

        final_set = []
        for f in range(len(seq_train_n)):
            final_set.append(list(np.hstack((seq_train_n[f],lsegs_n[f],vsegs_n[f],jsegs_n[f]))))
        final_set = np.array(final_set)

    return (final_set, posl)

def translateV(gene):
    if 'v' in gene.lower():
        v_gene = gene.split('*')[0]     
        return gene_to_num_str(v_gene,'V')
    else:
        return gene
    
def translateJ(gene):
    if 'j' in gene.lower():
        j_gene = gene.split('*')[0]     
        return gene_to_num_str(j_gene,'J')
    else:
        return gene
    
def produce_freqs(full_table):
    aa = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V',  'W', 'Y']
    nn=3 ## checked that stable wrt to nn in [2,3,4]
    fragments=[]
    for t in range(len(full_table)):
        lr = int(len(full_table[t][0])/2)
        ll= len(full_table[t][0]) - int(len(full_table[0][0])/2)
        fragments += (full_table[t][0][lr-nn:lr] + full_table[t][0][ll:ll+nn])
    freqs=[]
    for y in range(len(aa)):
        freqs.append(fragments.count(aa[y])/len(fragments))

    return freqs  

def convert_letterLR(final_set, SA, ltypes, posl, freqs):
    
    import random
    
    count_right=0 
    count_less=0 
    count_more=0
    ind_right=[]
    ind_less=[] 
    ind_more=[]
    aa = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V',  'W', 'Y']
    aadictinv = {k: aa[k] for k in range(len(aa))}
    final_set_converted=[]
    for u in range(len(final_set)):
        seq = final_set[u][:posl]
        lenn = final_set[u][posl]
        lennF = int(ltypes[lenn][1:])
        
        seqt = []
        for i in range(SA):
            if seq[i] != len(aa):
                seqt += aadictinv[seq[i]]

        seqt_m = []
        for i in range(-SA,0):
            if seq[i] != len(aa):
                seqt_m += aadictinv[seq[i]]
                
        seqt_fin = ''.join(seqt + seqt_m)
        
        if len(seqt_fin) == lennF:
            final_set_converted.append(seqt_fin)
            count_right+=100
            ind_right.append(u)
        
        if len(seqt_fin) < lennF:
            lendiff = lennF - len(''.join(seqt) + ''.join(seqt_m))
            #add_sym = random.sample(list(np.arange(len(aa))),lendiff) ## Here I should sample them with frequency of aa in CDR3
            rng = np.random.default_rng()
            add_sym = rng.choice(len(aa), lendiff, p = freqs)  
            final_set_converted.append(''.join(seqt) + ''.join([aadictinv[i] for i in add_sym]) + ''.join(seqt_m))
            count_less+=1
            ind_less.append(u)
        if len(seqt_fin) > lennF:
            lendiff = len(''.join(seqt) + ''.join(seqt_m)) - lennF
            count_more+=1
            ind_more.append(u)
            diff = int(lendiff/2)
            p = random.uniform(0, 1)
            if p<=0.5: ## checked that deleting on left/right or randomly alternating is completely equivalent at the level of pgen distribution
                if diff != 0:
                    final_set_converted.append(''.join(seqt[:-diff]) + ''.join(seqt_m[(lendiff - diff):]))
                else:
                    final_set_converted.append(''.join(seqt) + ''.join(seqt_m[(lendiff - diff):]))
            else:
                final_set_converted.append(''.join(seqt[:-(lendiff - diff)]) + ''.join(seqt_m[diff:]))
        
    return (final_set_converted, count_right, count_less, count_more, ind_right, ind_less, ind_more)
    



if makeLR:
    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
else:
    amino_acids = 'ACDEFGHIKLMNPQRSTVWY-'


if makeVJ:
    strVJ = 'withVJ'
    
else:
    strVJ = 'withoutVJ'
    
add_vj=strVJ 

add_lr=''
if makeLR:
    add_lr = '_LR'
str_lr = add_lr

RR=1

SA=20
SA1=20
if alpha_chain:
   SA=19
   SA1=19
if makeLR:
    SA=8

add_alpha = ''
chain_name = 'B'
if alpha_chain:
    add_alpha = '_alpha'
    chain_name='A'

if not(makeVJ):
    vtypes=[]
    jtypes=[]
        
if not(makeLR):
    ltypes=[]
    posl=-1

dataset = 'diffRBM_'+ pep ## This is the name of the folder with the model; if does not exist I create it 
name_folder = rootf + '/Example_diffRBM/'

t_df = pd.read_csv(name_folder + 'J2num_dict_train_train' + add_alpha + '.csv')
listj = list(t_df['J_gene'].values)
J2num_dict = dict(zip(t_df['J_gene'], t_df['RBM_VU_id']))
t_df = pd.read_csv(name_folder + 'V2num_dict_train_train' + add_alpha + '.csv')
listv = list(t_df['V_gene'].values)
V2num_dict = dict(zip(t_df['V_gene'], t_df['RBM_VU_id']))

vdgdb_df = pd.read_csv(name_folder + 'training_data/' + pep + '_aligned_' +str(SA1) +  add_alpha + '.csv')
Vlist = vdgdb_df['TR' + chain_name + 'V_gene'].to_list()
Jlist = vdgdb_df['TR' + chain_name + 'J_gene'].to_list()
inds = [t for t in range(len(vdgdb_df)) if (vdgdb_df['TR' + chain_name + 'V_gene'].values)[t] in listv and (vdgdb_df['TR' + chain_name + 'J_gene'].values)[t] in listj]
vdgdb_df = vdgdb_df.iloc[inds]
vdgdb_df = vdgdb_df[vdgdb_df['Label']==1]
train_data = vdgdb_df.drop_duplicates().reset_index(drop=True)
t_trainseq = train_data["ali_seq"].to_list()
if alpha_chain==0:
    seqs = train_data['CDR3_beta'].to_list()
if alpha_chain:
    seqs = train_data['CDR3_alpha'].to_list()
trainseq_2num = np.array(Proteins_utils.seq2num(t_trainseq), dtype=np.int16)
Vlist = train_data['TR' + chain_name + 'V_gene'].to_list()
Jlist = train_data['TR' + chain_name + 'J_gene'].to_list()
train_2num_withVJ = add_VJ_info_2num(trainseq_2num, Vlist, Jlist, V2num_dict, J2num_dict)
if makeVJ:
    train_dataI = np.copy(train_2num_withVJ)
else:
    train_dataI = np.copy(train_2num_withVJ[:,:-2])
if makeLR:
    full_table=[]
    full_tableS=[]
    for i in range(len(seqs)):
        if len(seqs[i])>=5 and len(seqs[i])<=30:
            full_table.append([seqs[i], t_trainseq[i], Vlist[i], Jlist[i]])
            full_tableS.append([seqs[i], Vlist[i], Jlist[i]])
    
    (train_dataI, posl) = encodeLR(full_tableS, SA, makeVJ)
    
    def convert_letterL(seqs_n):
        ldictinv = {k: ltypes[k] for k in range(len(ltypes))} 
        seqs = [ldictinv[e] for e in seqs_n]
        return seqs
    def convert_numberL(seqs): # convert to numbers the L types
        ldict = {ltypes[k]: k for k in range(len(ltypes))}
        l_num = np.array(list(map(lambda x: [ldict[y] for y in x], seqs[0:])), dtype=curr_int, order="c")
        return l_num

    if posl==-1 and makeVJ:
        print('check consistency VJ types')
        
        
vdgdb_df = pd.read_csv(name_folder + 'test_data/' + pep + '_aligned_' +str(SA1) +  add_alpha + '.csv')
Vlist = vdgdb_df['TR' + chain_name + 'V_gene'].to_list()
Jlist = vdgdb_df['TR' + chain_name + 'J_gene'].to_list()
inds = [t for t in range(len(vdgdb_df)) if (vdgdb_df['TR' + chain_name + 'V_gene'].values)[t] in listv and (vdgdb_df['TR' + chain_name + 'J_gene'].values)[t] in listj]
test_data_df = vdgdb_df.iloc[inds]
t_trainseq = test_data_df["ali_seq"].to_list()
if alpha_chain==0:
    seqs = test_data_df['CDR3_beta'].to_list()
if alpha_chain:
    seqs = test_data_df['CDR3_alpha'].to_list()
trainseq_2num = np.array(Proteins_utils.seq2num(t_trainseq), dtype=np.int16)
Vlist = test_data_df['TR' + chain_name + 'V_gene'].to_list()
Jlist = test_data_df['TR' + chain_name + 'J_gene'].to_list()
train_2num_withVJ = add_VJ_info_2num(trainseq_2num, Vlist, Jlist, V2num_dict, J2num_dict)
if makeVJ:
    test_data = np.copy(train_2num_withVJ)
else:
    test_data = np.copy(train_2num_withVJ[:,:-2])
if makeLR:
    full_table=[]
    full_tableS=[]
    for i in range(len(seqs)):
        if len(seqs[i])>=5 and len(seqs[i])<=30:
            full_table.append([seqs[i], t_trainseq[i], Vlist[i], Jlist[i]])
            full_tableS.append([seqs[i], Vlist[i], Jlist[i]])
    
    (test_data, posl_test) = encodeLR(full_tableS, SA, makeVJ)

    if posl_test==-1 and makeVJ:
        print('check consistency VJ type')

if os.path.exists(name_folder) is False:
    os.mkdir(name_folder)
if os.path.exists(name_folder + dataset) is False:
    os.mkdir(name_folder+dataset)

name_mf = '/models'
if os.path.exists(name_folder + dataset+ name_mf + '/') is False:
    os.mkdir(name_folder + dataset+ name_mf + '/')


if alpha_chain==0:
    Mb=1000000
    filename_cdr3 = name_folder + 'Background_models/train_data_0_aligned_20.txt' 
    filename_cdr3raw = name_folder + 'Background_models/train_data_0.txt' 
    train_seq0=[]
    with open(filename_cdr3) as f:
        for line in f:
            linesplit = line.strip().split('\n')
            train_seq0.append(linesplit[0])
    train_dataB_num = Proteins_utils.seq2num(train_seq0)                       
    if makeVJ:
        
        ffraw = pd.read_csv(filename_cdr3raw,sep='\t',header=None)
        seqs=list(np.array(ffraw)[:,0])
        VlistB = list(np.array(ffraw)[:,1])
        JlistB = list(np.array(ffraw)[:,2])
        train_2num_withVJ_B = add_VJ_info_2num(train_dataB_num, VlistB, JlistB, V2num_dict, J2num_dict)
        train_dataBI = train_2num_withVJ_B[:Mb]
    else:
        train_dataBI = train_dataB_num[:Mb]
       
if alpha_chain==1:
    df=pd.read_csv(name_folder + 'Background_models/alpha_background_aligned_19.csv')
    list_ali = list(df['ali_seq'].values)
    ind_list = [p for p in range(len(list_ali)) if list_ali[p][-1]!='-']
    df = df.iloc[ind_list]
    train_seq0 = list(df['ali_seq'])
    seqs = list(df['CDR3_alpha'])
    train_dataB_num = Proteins_utils.seq2num(train_seq0) 
    Mb = len(df)
    if makeVJ:  
        ffraw=df
        VlistB = list(np.array(ffraw)[:,0])
        JlistB = list(np.array(ffraw)[:,2])
        train_2num_withVJ_B = add_VJ_info_2num(train_dataB_num, VlistB, JlistB, V2num_dict, J2num_dict)
        train_dataBI = train_2num_withVJ_B
    else:
        train_dataBI = train_dataB_num

    
if makeLR:
    full_table=[]
    full_tableS=[]
    for i in range(len(seqs)):
        if len(seqs[i])>=5 and len(seqs[i])<=30:
            full_table.append([seqs[i], train_seq0[i], VlistB[i], JlistB[i]])
            full_tableS.append([seqs[i], VlistB[i], JlistB[i]])
    
    (train_dataBI, poslB) = encodeLR(full_tableS, SA, makeVJ)
    

if makeVJ:
    n_cv1 = np.max(np.unique(train_dataBI[:,-2])) + 1
    n_cv2 = np.max(np.unique(train_dataBI[:,-1])) + 1
    n_cv=np.max([n_cv1,n_cv2])
else:
    n_cv=21
    

listregB = [0.001] ## lists of weight regularization in background RBM
listhuB = [100]## lists of hidden units in background RBM

## Various parameters ##
th = 0.15 ## similarity threshold if you set the reweighting
## RBM parameters ##
ZF = False ## control the introduction of fields ##
BN = True ## control batch_size ##

visible = 'Potts' # Nature of visible units potential. Here, Potts states.
hidden = 'dReLU' # Nature of hidden units potential. Here, dReLU potential.
seed = 0
decay_after = 0.5 # Decay learning rate after 50% of iterations (default: 0.5). Value for RBM shown in paper: 0.5
N_MC = 15



## percentage to use in the training dataset ##
B = 100 ## percentage for positives ##
BB = 100 ## percentage for background ##

for lib in listreg:
    for hh in listhu:
        for hhB in listhuB:
            for libB in listregB:
                
                for repl in range(RR):

                    
                    seqs_temp = list(train_dataI)
                    
                    full_int = list(np.arange(int(len(seqs_temp)*1)))
                    full_intR = random.sample(full_int, int(1*(len(seqs_temp)))) 
                    data = np.array([seqs_temp[t] for t in full_intR])


                    train_data = data[:int(B*len(data)/100)]

                    if makeRW:
                        reweight = utilities_diffrbm.compute_MSA_weights(train_data, threshold = th)
                     

                    Mt = len(train_data)

                    train_dataB = train_dataBI[:int(BB*len(train_dataBI)/100)]
                   
                    Mb = len(train_dataB)

                    if makeRW:
                        reweightB = utilities_diffrbm.compute_MSA_weights(train_dataB, threshold = th)
                    
                    reg_fields_diffpos = 1/Mt
                    reg_weights_diffpos = lib

                    n_v = train_dataB.shape[1] # Number of visible units = # sites in alignment.

                    l2f = 1/Mb

                    n_hB = hhB
                    l1bB = libB
                    decay_after = 0.5
                    N_MC = 10

                    batch_size = 2000
                    n_iter = 40

                    out_par0 = '_RW' + str(makeRW) + '_TR' + str(BB) + add_vj + '_' + str(Mb) + add_alpha + str_lr
                    out_par = '_nh' + str(n_hB) + '_l12' + str(libB) + '_ZF' + str(ZF) + out_par0
                    nameB = name_folder + 'Background_models/backRBM' + out_par +'.data'
                    

                    if maketrainingB:
                        RBM_back = rbm.RBM(n_h = hhB, n_v = len(train_dataB[0]),n_cv= n_cv, visible='Potts', hidden='dReLU',random_state = seed, zero_field = False)
                        if maketrainingB:
                            RBM_back.fit(train_dataB, weights = None, batch_size = batch_size, n_iter = n_iter, l1b = libB, l2_fields = l2f, N_MC = N_MC, decay_after = decay_after, verbose = 1, shuffle_data=False, CD=False)
                            RBM_utils.saveRBM(nameB, RBM_back)

                    RBM_back =  RBM_utils.loadRBM(nameB)
                    

                    l2f = reg_fields_diffpos
                    

                    learning_rate = None # default behaviour
                    decay_after = 0.5
                    N_MC = 10
                    
                    if len(train_data) < 500:
                        bb = 16
                    else:
                        bb = 100
                    if len(train_data) > 3000: 
                        bb = 200

                    nni = int(2e4) // (train_data.shape[0] // bb)
                    nni=1
                    
                    batch_size = bb
                    n_iter = nni

                    n_v = train_data.shape[1] # Number of visible units = # sites in alignment.


                    
                    if hh == 0:
                        name_top = name_folder + dataset+ name_mf + '/topmodellin_imm'
                        
                        n_h_top = hh + 10

                        l2f = reg_fields_diffpos

                        batch_size = bb
                        n_iter = nni
                        
                        out_par0 = '_RW' + str(makeRW) + '_TR' + str(B) + add_vj + add_alpha + str_lr
                        out_par = '_lin_ZF' + str(ZF) + out_par0

                        RBMpost_lin = rbm.RBM(n_h = RBM_back.n_h + n_h_top, n_v = RBM_back.n_v, n_cv= RBM_back.n_cv, visible=visible, hidden=hidden, zero_field=ZF)
                        dRBM_lin = diffrbm.DiffRBM(RBM_back, RBMpost_lin)
                        dRBM_lin.update_post_from_back(vlayer=True, hlayer=True)

                        if maketraining:
                            if makeRW:
                                dRBM_lin.fit_top(train_data, weights=reweight, n_iter=n_iter, batch_size=batch_size, l2_fields = l2f, l1b=reg_weights_diffpos, N_MC=N_MC, decay_after=decay_after, verbose=vv, vverbose=vv, batch_norm=BN, callback = produce_callback_weights(dRBM_lin))
                            else:
                                dRBM_lin.fit_top(train_data, weights=None, n_iter=n_iter, batch_size=batch_size, l2_fields = l2f, l1b=reg_weights_diffpos, N_MC=N_MC, decay_after=decay_after, verbose=vv, vverbose=vv, batch_norm=BN, callback = produce_callback_weights(dRBM_lin))
                            RBM_utils.saveRBM(name_top + out_par + '.data', dRBM_lin)
                        dRBM_lin =  RBM_utils.loadRBM(name_top + out_par + '.data')
                        
                        ll = -dRBM_lin.top_rbm().free_energy(test_data.astype(np.int16))

                    else:
                        name_top = name_folder + dataset+ name_mf +  '/topmodel_imm'
                        n_h_top = hh
                        
                        l2f = reg_fields_diffpos

                        batch_size = bb
                        n_iter = nni
                        
                        out_par0 = '_RW' + str(makeRW) + '_TR' + str(B) + add_vj + add_alpha + str_lr
                        out_par = '_nh' + str(n_h_top) + '_l12' + str(lib) + '_ZF' + str(ZF) + out_par0
                        
                        RBMpost = rbm.RBM(n_h = RBM_back.n_h + n_h_top, n_v = RBM_back.n_v, n_cv= RBM_back.n_cv, visible='Potts', hidden=hidden, zero_field=ZF)
                        dRBM = diffrbm.DiffRBM(RBM_back, RBMpost)
                        dRBM.update_post_from_back(vlayer=True, hlayer=True)
                       
                        if maketraining:
                            if makeRW:
                                dRBM.fit_top(train_data, weights=reweight, n_iter=n_iter, batch_size=batch_size, l2_fields = l2f, l1b=reg_weights_diffpos, N_MC=N_MC, decay_after=decay_after, verbose=vv, vverbose=vv, batch_norm=BN)
                            else:
                                dRBM.fit_top(train_data, weights=None, n_iter=n_iter, batch_size=batch_size, l2_fields = l2f, l1b=reg_weights_diffpos, N_MC=N_MC, decay_after=decay_after, verbose=vv, vverbose=vv, batch_norm=BN)
                            RBM_utils.saveRBM(name_top + out_par + '.data', dRBM)
                        dRBM =  RBM_utils.loadRBM(name_top + out_par + '.data')

                        
                        ll = -dRBM.top_rbm().free_energy(test_data.astype(np.int16))
                        
                    test_data_df['Scores'] = ll
                    test_data_df.to_csv(name_folder + dataset + '/test_data_scored' + out_par + '.csv', sep='\t')
