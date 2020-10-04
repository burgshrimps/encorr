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
            num_spikes_neuron = []
            for i in range(10):
                stoi = cut_out(neuron, P.cut_time_before_stim, P.cut_time_after_stim, P.ts_rndm_baseline[i])
                stoi_spktimes_neuron.append(stoi)
                num_spikes_neuron.append(len(stoi))
            stoi_spktimes_tet.append(stoi_spktimes_neuron)
            num_spikes_tet.append(num_spikes_neuron)

    if phase == 'study':
        for neuron in tet:
            stoi_spktimes_neuron = []
            num_spikes_neuron = []
            for i in range(10):
                stoi = cut_out(neuron, P.cut_time_before_stim, P.cut_time_after_stim, P.ts_stim_study[i])
                stoi_spktimes_neuron.append(stoi)
                num_spikes_neuron.append(len(stoi))
            stoi_spktimes_tet.append(stoi_spktimes_neuron)
            num_spikes_tet.append(num_spikes_neuron)

    elif phase == 'exp_old':
        for neuron in tet:
            stoi_spktimes_neuron = []
            num_spikes_neuron = []
            for i in range(10):
                stoi = cut_out(neuron, P.cut_time_before_stim, P.cut_time_after_stim, P.ts_stim_exp_old[i])
                stoi_spktimes_neuron.append(stoi)
                num_spikes_neuron.append(len(stoi))
            stoi_spktimes_tet.append(stoi_spktimes_neuron)
            num_spikes_tet.append(num_spikes_neuron)

    elif phase == 'exp_new':
        for neuron in tet:
            stoi_spktimes_neuron = []
            num_spikes_neuron = []
            for i in range(10):
                stoi = cut_out(neuron, P.cut_time_before_stim, P.cut_time_after_stim, P.ts_stim_exp_new[i])
                stoi_spktimes_neuron.append(stoi)
                num_spikes_neuron.append(len(stoi))
            stoi_spktimes_tet.append(stoi_spktimes_neuron)
            num_spikes_tet.append(num_spikes_neuron)

    return stoi_spktimes_tet, num_spikes_tet


def get_mean_spikecount_exp(spikecount_exp_old, spikecount_exp_new):
    """ Receives spikecounts from exp old and exp new phases and retuns the average of the two. """
    mean_spkcount_exp = []
    for neur_idx in range(len(spikecount_exp_old)):
        mean_spkcount_exp.append(list(np.array((np.array((spikecount_exp_old[neur_idx]) + np.array(spikecount_exp_new[neur_idx])) / 2), dtype=int)))
    return mean_spkcount_exp


def get_stoi_baseline_spkcount(tet, P, spikecount_exp_old, spikecount_exp_new):
    """ Cuts out spiketrains of interest for baseline period based on spikecount in exp phase. """
    mean_spkcount_exp = get_mean_spikecount_exp(spikecount_exp_old, spikecount_exp_new)
    baseline_end_time = P.ts_stim_study[0] - 500000
    stoi_spktimes_tet = []
    for i, neuron in enumerate(tet):
        stoi_spktimes_neuron = []
        if neuron[0] < baseline_end_time:
            baseline_end_idx = np.where(neuron < baseline_end_time)[0][-1]
            neuron_baseline = neuron[:baseline_end_idx+1]
            start_idx = 0
            for j in range(10):
                stoi = neuron_baseline[start_idx:start_idx + mean_spkcount_exp[i][j]]
                stoi_spktimes_neuron.append(stoi)
                start_idx = start_idx + mean_spkcount_exp[i][j]
            stoi_spktimes_tet.append(stoi_spktimes_neuron)
        else:
            stoi_spktimes_tet.append([[]] * 10)
    return stoi_spktimes_tet

