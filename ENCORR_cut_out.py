import numpy as np


def cut_out(spiketrain, pre_start_time, post_end_time, stim_ts):
    """ Cuts out STOI (spike train of interest) from input spiketrain such that STOI = spiketrain[stim_ts-pre_start_time:stim_ts+post_end_time] """
    if spiketrain[-1] >= stim_ts - pre_start_time:
        if spiketrain[0] <= stim_ts - pre_start_time:
            start_idx = np.where(spiketrain >= stim_ts - pre_start_time)[0][0]
            end_idx = np.where(spiketrain <= stim_ts + post_end_time)[0][-1]
        else:
            if spiketrain[0] <= stim_ts + post_end_time:
                start_idx = 0
                end_idx = np.where(spiketrain <= stim_ts + post_end_time)[0][-1]
            else:
                return []
    else:
        return []
    return spiketrain[start_idx:end_idx+1]


def get_stoi(tet, P, phase):
    """ Iterates over all neurons in one tetrode and cuts out STOIs. """
    stoi_spktimes_tet = []
    num_spikes_tet = []
    if phase == 'baseline':
        for neuron in tet:
            stoi_spktimes_neuron = []
            for i in range(10):
                stoi_spktimes_neuron.append(cut_out(neuron, P.cut_time_before_stim, P.cut_time_after_stim, P.ts_rndm_baseline[i]))
            stoi_spktimes_tet.append(stoi_spktimes_neuron)

    if phase == 'study':
        for neuron in tet:
            stoi_spktimes_neuron = []
            for i in range(10):
                stoi_spktimes_neuron.append(cut_out(neuron, P.cut_time_before_stim, P.cut_time_after_stim, P.ts_stim_study[i]))
            stoi_spktimes_tet.append(stoi_spktimes_neuron)

    elif phase == 'exp_old':
        for neuron in tet:
            stoi_spktimes_neuron = []
            num_spikes_neuron = []
            for i in range(10):
                stoi_spktimes_neuron.append(cut_out(neuron, P.cut_time_before_stim, P.cut_time_after_stim, P.ts_stim_exp_old[i]))
            stoi_spktimes_tet.append(stoi_spktimes_neuron)

    elif phase == 'exp_new':
        for neuron in tet:
            stoi_spktimes_neuron = []
            num_spikes_neuron = 0
            for i in range(10):
                stoi = cut_out(neuron, P.cut_time_before_stim, P.cut_time_after_stim, P.ts_stim_exp_new[i])
                stoi_spktimes_neuron.append(stoi)
                num_spikes_neuron + len(stoi)
            stoi_spktimes_tet.append(stoi_spktimes_neuron)
            print(num_spikes_neuron)

    return stoi_spktimes_tet

