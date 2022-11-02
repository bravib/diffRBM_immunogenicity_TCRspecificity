import matplotlib.pyplot as plt

rootf = '.' ## Your folder ##

import sys
sys.path.append(rootf + '/diffRBM/source/')
sys.path.append(rootf + '/diffRBM/utilities/')

import argparse
parser = argparse.ArgumentParser() 
parser.add_argument('-pep', type = str, required=True)
parser.add_argument('-alpha', type = int, default = 0, required=False, help='Choose CDR3 on alpha chain instead of beta chain')
args = parser.parse_args()

alpha_chain = 0

if args.alpha:
    alpha_chain = 1
  
pep=args.pep

## Functions ##

import numpy as np
import pandas as pd
import subprocess, os

curr_int = np.int16

SA=20
if alpha_chain:
   SA=19

add_alpha = ''
chain_name = 'B'
if alpha_chain:
    add_alpha = '_alpha'
    chain_name='A'

name_folder = rootf + '/Example_diffRBM/'

if alpha_chain == 0:
    name_mat = name_folder + 'Align_utils/align_prot_to_seedhmmpy.py'
    SAmin = 3
    SAmax = 23

    filename_cdr3 = name_folder +'Background_models/train_data_0_aligned_20_ref.txt'
    seed=[]
    with open(filename_cdr3) as f:
        for line in f:
            linesplit = line.strip().split('\t')
            seed.append(linesplit[0])
            
    for task in ['training_data', 'test_data']:
    
        name = name_folder +task+'/' + pep +'.txt'

        selected = pd.read_csv(name, sep='\t')
        seqs_gapless = list(selected['TRB_CDR3'].values)

        Vseg = list(selected['TRBV'].values)
        Jseg = list(selected['TRBJ'].values)
        tcr_b =  list(seqs_gapless)

        if task != 'test_data':
            labels = list(selected['Label'].values)
        else:
            labels= np.ones((len(selected)))

        val_seqs=[]
        for tt in range(len(tcr_b)):
            if type(Vseg[tt]) != float and type(Jseg[tt]) != float:
                if '*' in Vseg[tt]:
                    vseg = Vseg[tt][:-3]
                if '*' in Jseg[tt]:
                    jseg = Jseg[tt][:-3]
                val_seqs.append([vseg,tcr_b[tt],jseg, labels[tt]])

        ## Align test data ##
        all_seqs1 = list(seed)
        name_seed = name_folder + 'Align_utils/prots_seed.txt'
        with open(name_seed, 'w') as out_f:
            for u in range(len(all_seqs1)):
                 
                out_f.write(all_seqs1[u] + '\n')

        name_seqs = name_folder + 'Align_utils/prots_seqs.txt'
        with open(name_seqs, 'w') as out_f:
            for u in range(len(seqs_gapless)):
                  
                out_f.write(seqs_gapless[u] + '\n')

        subprocess.call('python3 ' + name_mat + ' -sseed ' + name_seed + ' -sseqs ' + name_seqs + ' -SA ' + str(SA) + ' -SAmin ' + str(SAmin) + ' -SAmax ' + str(SAmax) + ' -yw 0', shell = True)
        seqs_al = np.loadtxt(name_folder + 'Align_utils/aligned_prot.txt')

        import os
        os.system('rm ' + name_folder + 'Align_utils/aligned_prot.txt')

        curr_int=np.int16
        def convert_letter(seqs_n): # convert to numbers already aligned seqs
            aa = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V',  'W', 'Y','-']
            aadictinv = {k: aa[k] for k in range(len(aa))} 
            seqs=[]
            if type(seqs_n[0]) == curr_int:
                seqs.append(''.join([aadictinv[e] for e in seqs_n]))
            else:
                for t in range(len(seqs_n)):
                    seqs.append(''.join([aadictinv[e] for e in seqs_n[t]]))
            return seqs

        sn = seqs_al.astype(np.int)
        seqs_g = convert_letter(sn)

        dataf = {'TRBV_gene': np.array(val_seqs)[:,0],
            'CDR3_beta': np.array(val_seqs)[:,1],
            'TRBJ_gene': np.array(val_seqs)[:,2],
            'Label': np.array(val_seqs)[:,3],
            'ali_seq': seqs_g
            }

        df = pd.DataFrame(dataf, columns = ['Label','TRBV_gene','CDR3_beta','TRBJ_gene','ali_seq'])
        df.to_csv(name_folder + task + '/' + pep + '_aligned_' + str(SA) + '.csv', index=False)
        
else:
    name_mat = name_folder + 'Align_utils/align_prot_to_seedhmmpy.py'
    SAmin = 5
    SAmax = 23
    
    for task in ['training_data', 'test_data']:
        name = name_folder +task+'/' + pep +'.txt'
       

        file=pd.read_csv(name_folder +'Background_models/alpha_background_aligned_' + str(SA) + '.csv')
        seed = list(file['ali_seq_ref'].values)


        name = name_folder +task+'/' + pep +'.txt'
        selected = pd.read_csv(name, sep='\t')
        seqs_gapless = list(selected['TRA_CDR3'].values)
        
        Vseg = list(selected['TRAV'].values)
        Jseg = list(selected['TRAJ'].values)
        tcr_b =  list(seqs_gapless)
        if task != 'test_data':
            labels = list(selected['Label'].values)
        else:
            labels = np.ones((len(selected)))

        val_seqs=[]
        for tt in range(len(tcr_b)):
            if type(Vseg[tt]) != float and type(Jseg[tt]) != float:
                vseg = Vseg[tt]
                if '*' in Vseg[tt]:
                    vseg = Vseg[tt][:-3]
                if '/' in vseg:
                    i = vseg.index('/')
                    vseg = vseg[:i]
                jseg = Jseg[tt]
                if '*' in Jseg[tt]:
                    jseg = Jseg[tt][:-3]
                if '/' in jseg:
                    i = jseg.index('/')
                    jseg = jseg[:i]
                val_seqs.append([vseg,tcr_b[tt],jseg,labels[tt]])
            


        ## Align test data ##
        all_seqs1 = list(seed)
        name_seed = name_folder + 'Align_utils/prots_seed.txt'
        with open(name_seed, 'w') as out_f:
            for u in range(len(all_seqs1)):
                  
                out_f.write(all_seqs1[u] + '\n')

        name_seqs = name_folder + 'Align_utils/prots_seqs.txt'
        with open(name_seqs, 'w') as out_f:
            for u in range(len(seqs_gapless)):
                  
                out_f.write(seqs_gapless[u] + '\n')

        subprocess.call('python3 ' + name_mat + ' -sseed ' + name_seed + ' -sseqs ' + name_seqs + ' -SA ' + str(SA) + ' -SAmin ' + str(SAmin) + ' -SAmax ' + str(SAmax) + ' -yw 0', shell = True)
        seqs_al = np.loadtxt(name_folder + 'Align_utils/aligned_prot.txt')

        import os
        os.system('rm ' + name_folder + 'Align_utils/aligned_prot.txt')

        curr_int=np.int16
        def convert_letter(seqs_n): # convert to numbers already aligned seqs
            aa = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V',  'W', 'Y','-']
            aadictinv = {k: aa[k] for k in range(len(aa))} 
            seqs=[]
            if type(seqs_n[0]) == curr_int:
                seqs.append(''.join([aadictinv[e] for e in seqs_n]))
            else:
                for t in range(len(seqs_n)):
                    seqs.append(''.join([aadictinv[e] for e in seqs_n[t]]))
            return seqs

        sn = seqs_al.astype(np.int)
        seqs_al = convert_letter(sn)
        
        inds = [s for s in range(len(seqs_al)) if seqs_al[s][-1]=='-']
        seqs_al_temp = list(seqs_al)
        for i in inds:
            j=seqs_al_temp[i].index('-')
            ss=list(seqs_al_temp[i])
            seqs_cut = seqs_gapless[i][j:]
            for f in range(1,len(seqs_cut)):
                ss[-f] = seqs_cut[-f]
            seqs_al_temp[i] = ''.join(ss)
        seqs_g=seqs_al_temp
        
        dataf = {'TRAV_gene': np.array(val_seqs)[:,0],
            'CDR3_alpha': np.array(val_seqs)[:,1],
            'TRAJ_gene': np.array(val_seqs)[:,2],
            'Label': np.array(val_seqs)[:,3],
            'ali_seq': seqs_g
            }
        
        df = pd.DataFrame(dataf, columns = ['Label','TRAV_gene','CDR3_alpha','TRAJ_gene','ali_seq'])
        df.to_csv(name_folder+task+'/' + pep + '_aligned_' + str(SA) + '_alpha.csv', index=False)

