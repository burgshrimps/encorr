import os
import numpy as np
import fnmatch as fnm
from scipy.io import loadmat

class Parameters:
    # Conatins all parameters of the experiment
    def __init__(self, tetrodes, labels, stim_study, stim_exp, stim_exp_new, stim_exp_old, ts_resp_exp, ts_resp_exp_new, ts_resp_exp_old,
            ts_stim_all, ts_stim_study, ts_stim_exp_old, ts_stim_exp_new, ts_stim_exp, cut_time_before_stim, cut_time_study_after_stim,
            cut_time_exp_after_resp):
        self.tetrodes = tetrodes  # list of brain regions associated whith each tetrode
        self.labels = labels  # labels for old/new odor during retrieval phase (1 = new, 0 = old)
        self.stim_study = stim_study  # number code for odors used in study phase
        self.stim_exp = stim_exp  # number code for odors used in retrieval phase
        self.stim_exp_new = stim_exp_new  # number code for new odors used in retrieval phase
        self.stim_exp_old = stim_exp_old  # number code for old odors used in retrieval phase
        self.ts_resp_exp = ts_resp_exp  # timestamps of rodent response to stimulus during retrieval phase in [ms]
        self.ts_resp_exp_new = ts_resp_exp_new  # timestamps of rodent response to new odors during retrieval phase in [ms]
        self.ts_resp_exp_old = ts_resp_exp_old  # timestamps of rodent response to old odors during retrieval phase in [ms]
        self.ts_stim_all = ts_stim_all  # timestamps of all stimulus presentations in [ms]
        self.ts_stim_study = ts_stim_study   # timestamps of stimulus presentation during study phase in [ms]
        self.ts_stim_exp_old = ts_stim_exp_old  # timestamps of stimulus presentation of old odors during retrieval phase in [ms]
        self.ts_stim_exp_new = ts_stim_exp_new  # timestamps of stimulus presentation of new odors during retrieval phase in [ms]
        self.ts_stim_exp = ts_stim_exp  # timestamps of stimulus presentation of all odors during retrieval phase in [ms]
        self.cut_time_before_stim = cut_time_before_stim  # time before stim. presentation timestamp = start of each spiketrain of interest (STOI)
        self.cut_time_study_after_stim = cut_time_study_after_stim # time after stim. presentation timestamp in study phase = end of each (STOI) in study phase
        self.cut_time_exp_after_resp = cut_time_exp_after_resp  # time after resp. timestamp = end of each STOI in exp phase


def loadparams(params_mat, cut_time_before_stim, cut_time_exp_after_resp):
    param = loadmat(params_mat)
    tetrodes = np.array([tet[0] for tet in param['tet_list'][0]]) 
    #labels = np.array([label[0] for label in param['label_oldnew']]) 
    labels = np.array(param['label_oldnew'][0]) ### changed for LE84
    stim_study = np.array([stim[0] for stim in param['stim_list_study']]) 
    stim_exp =  param['stim_list_exp'][0] 
    stim_exp_new = stim_exp[np.where(labels == 1)] 
    stim_exp_old = stim_exp[np.where(labels == 0)] 
    ts_resp_exp = param['ts_response_video'][:,0] ### chnaged for LE84, change back to [:,3] * 1000 for LE46
    ts_resp_exp_new = ts_resp_exp[np.where(labels == 1)] 
    ts_resp_exp_old = ts_resp_exp[np.where(labels == 0)] 
    ts_stim_all = [] ### changed for LE84, change back to param['ts_stimon'][0] for LE46
    ts_stim_study = param['ts_stimon_study'][:,0] ### changed for LE84, change back to ts_stim_all[0:10] for LE46
    ts_stim_exp = param['ts_stimon_exp'][:,0] ### changed for LE48, change back to param['ts_stimon_exp_video'][:,3] for LE46
    ts_stim_exp_new = ts_stim_exp[np.where(labels == 1)] 
    ts_stim_exp_old = ts_stim_exp[np.where(labels == 0)] 
    cut_time_study_after_stim = np.mean((ts_resp_exp+cut_time_exp_after_resp)-ts_stim_exp, dtype = int) # mean spiketrain of interest (STOI) length in exp phase
    return Parameters(tetrodes, labels, stim_study, stim_exp, stim_exp_new, stim_exp_old, ts_resp_exp, ts_resp_exp_new, ts_resp_exp_old,
                      ts_stim_all, ts_stim_study, ts_stim_exp_old, ts_stim_exp_new, ts_stim_exp, cut_time_before_stim, 
                      cut_time_study_after_stim, cut_time_exp_after_resp)


def loadtet(tet_mat_file, sampling_rate):
    tet_mat = loadmat(tet_mat_file)
    tet = [neuron[:,0].astype(int) for neuron in tet_mat['spikes']['spktimes'][0][0][0] / sampling_rate] 
    return tet