import numpy as np


def cut_out(spiketrain, pre_start_time, post_end_time, start_ts, end_ts):
    """ Cuts out STOI from input spiketrain such that STOI = spiketrain[start_ts-pre_start_time:end_ts+post_end_time] """
    if spiketrain[0] <= start_ts - pre_start_time:
        start_idx = np.where(spiketrain >= start_ts - pre_start_time)[0][0]
        end_idx = np.where(spiketrain <= end_ts + post_end_time)[0][-1]
    else:
        if spiketrain[0] <= end_ts + post_end_time:
            start_idx = 0
            end_idx = np.where(spiketrain <= end_ts + post_end_time)[0][-1]
        else:
            return []
    return spiketrain[start_idx:end_idx+1]


def get_stoi(tet, P, phase):
    """ Iterates over all neurons in one tetrode and cuts out STOIs. """
    stoi_spktimes_tet = []
    if phase == 'baseline':
        for neuron in tet:
            try:
                end_idx = np.where(neuron <= 300000)[0][-1]
            except IndexError:
                end_idx = 0
            stoi_spktimes_tet.append([neuron[:end_idx]])

    if phase == 'study':
        for neuron in tet:
            stoi_spktimes_neuron = []
            for i in range(10):
                stoi_spktimes_neuron.append(cut_out(neuron, P.cut_time_before_stim, P.cut_time_study_after_stim, P.ts_stim_study[i],
                                                    P.ts_stim_study[i]))
            stoi_spktimes_tet.append(stoi_spktimes_neuron)

    elif phase == 'exp_old':
        for neuron in tet:
            stoi_spktimes_neuron = []
            for i in range(10):
                stoi_spktimes_neuron.append(cut_out(neuron, P.cut_time_before_stim, P.cut_time_exp_after_resp, P.ts_stim_exp_old[i],
                                                    P.ts_resp_exp_old[i]))
            stoi_spktimes_tet.append(stoi_spktimes_neuron)

    elif phase == 'exp_new':
        for neuron in tet:
            stoi_spktimes_neuron = []
            for i in range(10):
                stoi_spktimes_neuron.append(cut_out(neuron, P.cut_time_before_stim, P.cut_time_exp_after_resp, P.ts_stim_exp_new[i],
                                                    P.ts_resp_exp_new[i]))
            stoi_spktimes_tet.append(stoi_spktimes_neuron)

    return stoi_spktimes_tet

