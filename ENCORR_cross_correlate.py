import elephant as eph
import quantities as pq
import neo
import numpy as np
import logging


class CorrelationRecord:

    def __init__(self, ref_tet, ref_neur, tar_tet, tar_neur, fmt, phases):
        self.ref_tet = ref_tet
        self.ref_neur = ref_neur
        self.tar_tet = tar_tet
        self.tar_neur = tar_neur
        self.format = fmt
        self.phases = phases


class CorrelationFile:

    def __init__(self, outfile, read_write):
        self.outfile = outfile
        self.f = open(self.outfile, read_write)

    def write(self, corr_rec):
        phases_as_string = ''
        for phase in corr_rec.phases:
            phases_as_string += np.array2string(phase['CH'], max_line_width=1000000) + ':' + str(phase['RS']) + '\t'
        corr_line = 'tet{0}\t{1}\ttet{2}\t{3}\t{4}\t{5}\n'.format(corr_rec.ref_tet, corr_rec.ref_neur+1, corr_rec.tar_tet, corr_rec.tar_neur+1, ':'.join(corr_rec.format), phases_as_string)
        self.f.write(corr_line)

    def parse_line(self, line):
        fields = line.split('\t')
        ref_tet = fields[0]
        ref_neur = fields[1]
        tar_tet = fields[2]
        tar_neur = fields[3]
        fmt = fields[4].split(':')
        phases_as_string = fields[5:]
        phases = []
        for phase_as_str in phases_as_string:
            phase_dict = dict()
            phase_dict['CH'] = np.fromstring(phase_as_str.split(':')[0][1:-1], sep=' ', dtype=int)
            phase_dict['RS'] = int(phase_as_str.split(':')[1])
            phases.append(phase_dict)
        return CorrelationRecord(ref_tet, ref_neur, tar_tet, tar_neur, fmt, phases)
        
    def fetch(self):
        for line in self.f.readlines():
            yield self.parse_line(line.strip())

            
def list2neo(st1, st2):
    """ Converts two lists of spike times into neo.Spiketrain objects. """
    t_stop = max(max(st1), max(st2)) + 1
    neo_st1 = neo.core.SpikeTrain(st1*pq.ms, t_stop=t_stop)
    neo_st2 = neo.core.SpikeTrain(st2*pq.ms, t_stop=t_stop)
    return neo_st1, neo_st2


def cch(st1, st2, binsize, windowsize, border_correction):
    """ Cross-correlates two neo spiketrains and returns cross-correlation histogram """
    st2.t_stop = st1.t_stop  # usually the two spiketrains cch() receives
                             # should have the same t_stop but sometimes a
                             # rounding error occured resulting in the two
                             # floating point numbers not being the same
    st1_binned = eph.conversion.BinnedSpikeTrain(st1, binsize=binsize * pq.ms)
    st2_binned = eph.conversion.BinnedSpikeTrain(st2, binsize=binsize * pq.ms)
    cch_object = eph.spike_train_correlation.cross_correlation_histogram(
        st1_binned, st2_binned, window=[-windowsize, windowsize],
        border_correction=border_correction, binary=False, kernel=None)
    cch = [int(val) for val in cch_object[0]]
    return np.array(cch, dtype=int)

def get_cch_for_all_neurons(ref_tet_id, tar_tet_id, ref_spiketimes, tar_spiketimes, phase, binsize, windowsize, border_correction):
    cch_all_neurons = dict()
    num_ref_spikes_all_neurons = dict()
    for rn in range(len(ref_spiketimes)):
        i = 1 if ref_tet_id == tar_tet_id else 0
        for tn in range(rn+i, len(tar_spiketimes)):
            logging.info('CCH @ {0}: RefNeur {1}/{2}, TarNeur {3}/{4}'.format(phase, rn+1, len(ref_spiketimes), tn+1, len(tar_spiketimes)))
            neuron_cch = np.zeros(2*windowsize+1, dtype=int)
            neuron_num_ref_spikes = 0
            for odor in range(len(ref_spiketimes[0])):
                try:
                    ref_spiketrain, tar_spiketrain = list2neo(ref_spiketimes[rn][odor], tar_spiketimes[tn][odor])
                    neuron_cch += cch(ref_spiketrain, tar_spiketrain, binsize, windowsize, border_correction)
                    neuron_num_ref_spikes += len(ref_spiketimes[rn][odor])
                except ValueError:  # case when at least one of the two spiketrains contains no neurons at all
                    neuron_cch += np.zeros(2*windowsize+1, dtype=int)
                    neuron_num_ref_spikes += len(ref_spiketimes[rn][odor])
            cch_all_neurons[(rn, tn)] = neuron_cch
            num_ref_spikes_all_neurons[(rn, tn)] = neuron_num_ref_spikes
    return cch_all_neurons, num_ref_spikes_all_neurons


def write_to_ccg(ref_tet_id, tar_tet_id, cch_baseline, cch_study, cch_exp_old, cch_exp_new, ref_spk_baseline, ref_spk_study,
                 ref_spk_exp_old, ref_spk_exp_new, outfile):
    ccg_out = CorrelationFile(outfile, 'w')
    fmt = ['CH', 'RS']
    for neur_pair in sorted(cch_baseline):
        phases = []

        baseline_info = dict()
        baseline_info['CH'] = cch_baseline[neur_pair]
        baseline_info['RS'] = ref_spk_baseline[neur_pair]
        phases.append(baseline_info)

        study_info = dict()
        study_info['CH'] = cch_study[neur_pair]
        study_info['RS'] = ref_spk_study[neur_pair]
        phases.append(study_info)

        exp_old_info = dict()
        exp_old_info['CH'] = cch_exp_old[neur_pair]
        exp_old_info['RS'] = ref_spk_exp_old[neur_pair]
        phases.append(exp_old_info)

        exp_new_info = dict()
        exp_new_info['CH'] = cch_exp_new[neur_pair]
        exp_new_info['RS'] = ref_spk_exp_new[neur_pair]
        phases.append(exp_new_info)

        corr_rec = CorrelationRecord(ref_tet_id, neur_pair[0], tar_tet_id, neur_pair[1], fmt, phases)
        ccg_out.write(corr_rec)

        

        