#import utilities as utilities
import rbm
import numpy as np
import pandas as pd
from numba import njit, prange
from scipy.optimize import linear_sum_assignment
from numpy import linalg as LA
import sys
sys.path.append('../utilities/')
import RBM_utils

curr_int = np.int16
curr_float = np.float32

def rbm_join_hidden(rbm1, rbm2):
    assert rbm1.n_v == rbm2.n_v
    assert rbm1.n_cv == rbm2.n_cv
    assert rbm1.n_ch == rbm2.n_ch
    assert rbm1.visible == rbm2.visible and rbm1.hidden == rbm2.hidden
    rbm_join = rbm.RBM(n_v=rbm1.n_v, n_h = rbm1.n_h + rbm2.n_h,
                       n_cv=rbm1.n_cv, n_ch=rbm1.n_ch,
                       visible=rbm1.visible, hidden=rbm1.hidden)
    
    if rbm_join.n_cv > 1 and rbm_join.n_ch > 1:
        rbm_join.weights[:rbm1.n_h,:,:,:] = rbm1.weights
        rbm_join.weights[rbm1.n_h:,:,:,:] = rbm2.weights
    elif rbm_join.n_cv > 1 or rbm_join.n_ch > 1:
        rbm_join.weights[:rbm1.n_h,:,:] = rbm1.weights
        rbm_join.weights[rbm1.n_h:,:,:] = rbm2.weights
    elif rbm_join.n_cv == rbm_join.n_ch == 1:
        rbm_join.weights[:rbm1.n_h,:] = rbm1.weights
        rbm_join.weights[rbm1.n_h:,:] = rbm2.weights
    else:
        raise Exception('Check n_cv, n_ch.')
    
    for key in rbm_join.vlayer.list_params:
        rbm_join.vlayer.__dict__[key] = rbm1.vlayer.__dict__[key] + rbm2.vlayer.__dict__[key]
    
    for key in rbm_join.hlayer.list_params:
        if rbm_join.n_ch == 1:
            rbm_join.hlayer.__dict__[key][:rbm1.n_h] = rbm1.hlayer.__dict__[key]
            rbm_join.hlayer.__dict__[key][rbm1.n_h:] = rbm2.hlayer.__dict__[key]
        else:
            rbm_join.hlayer.__dict__[key][:rbm1.n_h,:] = rbm1.hlayer.__dict__[key]
            rbm_join.hlayer.__dict__[key][rbm1.n_h:,:] = rbm2.hlayer.__dict__[key]

    rbm_join.moments_data = rbm1.moments_data
    return rbm_join


def rbm_extract_hidden(RBM, hid):
    sub_rbm = rbm.RBM(n_v=RBM.n_v, n_h=len(hid),
                      n_cv=RBM.n_cv, n_ch=RBM.n_ch,
                      visible=RBM.visible, hidden=RBM.hidden)
    
    if RBM.n_cv > 1 and RBM.n_ch > 1:
        sub_rbm.weights = RBM.weights[hid,:,:,:]
    elif RBM.n_cv > 1 or RBM.n_ch > 1:
        sub_rbm.weights = RBM.weights[hid,:,:]
    elif RBM.n_cv == rbm.n_ch == 1:
        sub_rbm.weights = RBM.weights[hid,:]
    else:
        raise Exception('Check n_cv, n_ch.')
    
    for key in RBM.vlayer.list_params:
        sub_rbm.vlayer.__dict__[key] = RBM.vlayer.__dict__[key]

    for key in RBM.hlayer.list_params:
        if RBM.n_ch == 1:
            sub_rbm.hlayer.__dict__[key] = RBM.hlayer.__dict__[key][hid]
        else:
            sub_rbm.hlayer.__dict__[key] = RBM.hlayer.__dict__[key][hid,:]
    
    return sub_rbm


def rbm_split_hidden(RBM, hid):
    hid_ = [u for u in range(RBM.n_h) if u not in hid]
    rbm1 = rbm_extract_hidden(RBM, hid)
    rbm2 = rbm_extract_hidden(RBM, hid_)
    return rbm1, rbm2


def weights_matching(rbm1, rbm2, match_overlaps=False):
    """
    Match the weights of rbm1 with those of rbm2, by maximizing the sum of the modulus
    of the scalar products of matched pairs.
    Return the matched pairs as a pandas dataframe with fields: number weight rbm1, 
    number weight rbm2, overlap, norm of weight 1, norm of weight 2; it also gives the
    together with the matching "cost" normalized so that it is between 0 and 1.
    If match_overlaps, (modulus of) overlaps are used instead of scalar products in 
    the cost matrix.
    """
    
    overlap_matrix = np.array([[np.dot(x.flatten(), y.flatten()) / (LA.norm(x)*LA.norm(y)) for x in rbm2.weights] for y in rbm1.weights])
    norms1 = np.array([LA.norm(x) for x in rbm1.weights])
    norms2 = np.array([LA.norm(x) for x in rbm2.weights])
    if not match_overlaps:
        cost_matrix = np.array([[np.abs(np.dot(x.flatten(), y.flatten())) for x in rbm2.weights] for y in rbm1.weights])
    else:
        cost_matrix = np.abs(overlap_matrix)
        
    row_ind, col_ind = linear_sum_assignment(cost_matrix, maximize=True)
    if not match_overlaps:
        total_cost = cost_matrix[row_ind, col_ind].sum() / (np.sum(norms1[row_ind] * norms2[col_ind]))
    else:
        total_cost = cost_matrix[row_ind, col_ind].sum() / len(cost_matrix[row_ind, col_ind])
    overlap_matching = overlap_matrix[row_ind, col_ind]
    
    matching_pairs_costs = pd.DataFrame({'rbm1_weight_num' : row_ind, 'rbm2_weight_num' : col_ind,
                                         'overlap' : overlap_matching, 'rbm1_weight_norm' : norms1[row_ind],
                                         'rbm2_weight_norm' : norms2[col_ind]})    
    return matching_pairs_costs, total_cost


def KL_divergence(rbm1, rbm2, beta=1, Nchains=2000, Nsequences=20000, Ntherm=1000, Nsteps=200):
    """
    Sample from rbm1 to estimate the KL divergence of the probability distribution for the visible
    units.
    """
    Lchains = int(Nsequences/Nchains) 
    sample_v, sample_h = RBM_utils.gen_data_lowT(rbm1, Nchains = Nchains, beta=beta, Lchains = Lchains, 
                                             Nstep=Nsteps, N_PT=1, update_betas=False, Nthermalize=Ntherm)
    return np.average(rbm1.likelihood((sample_v)) - rbm2.likelihood((sample_v))) 


@njit(parallel=False, cache=True, nogil=False)
def compute_MSA_weights(MSA, threshold = 0.1, verbose=False): 
    """
    Compute the MSA sequence weights as the inverse of the 
    number of neighbouring sequences.
    A neighbour of a sequence is another sequence within a given 
    threshold of Hamming distance (given as a fraction of the total 
    sequence length). MSA string sequences must be converted into 
    number sequences with seq2num function before passing it to 
    compute_MSA_weights.
    """
    
    B = MSA.shape[0]
    N = MSA.shape[1]
    num_neighbours = np.zeros(B)
    for a in range(0, B):
        num_neighbours[a] += 1
        t = MSA[a]
        if (a % 10000 == 0) and verbose:
            print(a)
        for b in prange(a + 1, B):
            dist = np.mean(t != MSA[b])
            if dist < threshold:
                num_neighbours[a] += 1
                num_neighbours[b] += 1
    return 1.0/num_neighbours

