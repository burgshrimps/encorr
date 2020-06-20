import os
import numpy as np
import fnmatch as fnm
from scipy.io import loadmat
import pandas as pd

class Parameters:
    """ Conatins all parameters of the experiment """
    def __init__(self, tetrodes, labels, ts_resp_exp, ts_resp_exp_new, ts_resp_exp_old, ts_stim_study, ts_stim_exp_old, ts_stim_exp_new, 
                 ts_stim_exp, cut_time_before_stim, cut_time_after_stim, ts_rndm_baseline):
        self.tetrodes = tetrodes  # list of brain regions associated whith each tetrode
        self.labels = labels  # labels for old/new odor during retrieval phase (1 = new, 0 = old)
        self.ts_resp_exp = ts_resp_exp  # timestamps of rodent response to stimulus during retrieval phase in [ms]
        self.ts_resp_exp_new = ts_resp_exp_new  # timestamps of rodent response to new odors during retrieval phase in [ms]
        self.ts_resp_exp_old = ts_resp_exp_old  # timestamps of rodent response to old odors during retrieval phase in [ms]
        self.ts_stim_study = ts_stim_study   # timestamps of stimulus presentation during study phase in [ms]
        self.ts_stim_exp_old = ts_stim_exp_old  # timestamps of stimulus presentation of old odors during retrieval phase in [ms]
        self.ts_stim_exp_new = ts_stim_exp_new  # timestamps of stimulus presentation of new odors during retrieval phase in [ms]
        self.ts_stim_exp = ts_stim_exp  # timestamps of stimulus presentation of all odors during retrieval phase in [ms]
        self.cut_time_before_stim = cut_time_before_stim  # time before stim. presentation timestamp = start of each spiketrain of interest (STOI)
        self.cut_time_after_stim = cut_time_after_stim # time after stim. presentation timestamp = end of each STOI
        self.ts_rndm_baseline = ts_rndm_baseline


def random_sample_baseline_ts(baseline_end_ts, cut_time_before_stim, cut_time_after_stim):
    """ Picks 10 random timepoints between 0 and the end of baseline period. The minimum interval between these
    timepoints is governed by cut_time_after_stim and cut_time_before_stim """
    baseline_ts = [np.random.randint(0,baseline_end_ts)]
    while len(baseline_ts) <= 10:
        rndm = np.random.randint(0,baseline_end_ts) 
        repeat = False
        for ts in baseline_ts:
            if abs(rndm - ts) < cut_time_after_stim + cut_time_before_stim:
                repeat = True
        if not repeat:
            baseline_ts.append(rndm)
    return sorted(baseline_ts)


def loadparams(params_mat, cut_time_before_stim, cut_time_after_stim, baseline_end_time):
    """ Loads parameters from .MAT file and creates instance of class Parameters """
    param = loadmat(params_mat)
    tetrodes = np.array([tet[0] for tet in param['tet_list'][0]]) 
    labels = np.array(param['label_oldnew'][0])
    ts_resp_exp = param['ts_response_video'][:,0] * 1000
    ts_resp_exp_new = ts_resp_exp[np.where(labels == 1)] 
    ts_resp_exp_old = ts_resp_exp[np.where(labels == 0)] 
    ts_stim_study = param['ts_stimon_study'][:,0] * 1000
    ts_stim_exp = param['ts_stimon_exp'][:,0] * 1000
    ts_stim_exp_new = ts_stim_exp[np.where(labels == 1)] 
    ts_stim_exp_old = ts_stim_exp[np.where(labels == 0)] 
    ts_rndm_baseline = random_sample_baseline_ts(baseline_end_time, cut_time_before_stim, cut_time_after_stim)

    return Parameters(tetrodes, labels, ts_resp_exp, ts_resp_exp_new, ts_resp_exp_old, ts_stim_study, ts_stim_exp_old, ts_stim_exp_new, 
                      ts_stim_exp, cut_time_before_stim, cut_time_after_stim, ts_rndm_baseline)


def loadtet(tet_mat_file, sampling_rate):
    """ Loads spike time information for a whole tetrode from .MAT file and normalises the 
    timestamps by dividing by the sampling rate """
    tet_mat = loadmat(tet_mat_file)
    tet = [np.rint(neuron[:,0]).astype(int) for neuron in tet_mat['spikes']['spktimes'][0][0][0] / (sampling_rate / 1000)] # first round to nearest int then convert to int 
    return tet

def load_tet_info(tet_info_mat):
    """ Loads tet_info .MAT file """
    tet_info = loadmat(tet_info_mat)
    areas = [area_array[0] for area_array in tet_info['tet_info'][0]]
    neur_counts = [neur_array[0][0] for neur_array in tet_info['tet_info'][1]]
    tet_info_df = pd.DataFrame(list(zip(areas, neur_counts)), columns=['area', 'neur_count'])
    return tet_info_df
