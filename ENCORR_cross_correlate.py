import elephant as eph
import quantities as pq
import neo
import numpy as np
import logging
import os


class CorrelationRecord:
    """ One record corresponds to one row in the CCG file. It contains the cross-correlation histograms
    between two neurons for all phases """

    def __init__(self, ref_area, ref_tet, ref_neur, tar_area, tar_tet, tar_neur, fmt, phases):
        self.ref_area = ref_area
        self.ref_tet = ref_tet
        self.ref_neur = ref_neur
        self.tar_area = tar_area
        self.tar_tet = tar_tet
        self.tar_neur = tar_neur
        self.format = fmt
        self.phases = phases


class CorrelationHeader:
    """ Header information for CCG file """

    def __init__(self, ref_mat, tar_mat, sampling_rate, cut_time_before_stim, cut_time_after_stim, binsize, windowsize, border_correction, fmt, col_names):
        self.ref_mat = ref_mat
        self.tar_mat = tar_mat
        self.sampling_rate = sampling_rate
        self.cut_time_before_stim = cut_time_before_stim
        self.cut_time_after_stim = cut_time_after_stim
        self.binsize = binsize
        self.windowsize = windowsize // binsize
        self.border_correction = border_correction
        self.format = fmt
        self.col_names = col_names


class CorrelationFile:
    """ CCG file class handling all reading and writing of CCG files """

    def __init__(self, file, read_write, header=None):
        self.header = header

        path_to_file = os.path.dirname(file)
        if not os.path.exists(path_to_file) and path_to_file != '':
            os.makedirs(path_to_file)
        self.f = open(file, read_write)

        if read_write == 'w' and self.header != None:
            self.write_header()
        
        elif read_write == 'r':
            self.all_lines = self.f.readlines()
            self.header = self.read_header()
            self.f.close()

    def read_header(self):
        fmt = []
        for line in self.all_lines:
            if line.startswith('#'):
                if line.startswith('##'):
                    line_splitted = line.strip().split('=')
                    if line_splitted[0][2:] == 'ref_mat':
                        ref_mat = line_splitted[1]
                    elif line_splitted[0][2:] == 'tar_mat':
                        tar_mat = line_splitted[1]
                    elif line_splitted[0][2:] == 'sampling_rate':
                        sampling_rate = float(line_splitted[1])
                    elif line_splitted[0][2:] == 'cut_time_before_stim':
                        cut_time_before_stim = int(line_splitted[1])
                    elif line_splitted[0][2:] == 'cut_time_after_stim':
                        cut_time_after_stim = int(line_splitted[1])
                    elif line_splitted[0][2:] == 'binsize':
                        binsize = int(line_splitted[1])
                    elif line_splitted[0][2:] == 'windowsize':
                        windowsize = int(line_splitted[1])
                        print(windowsize)
                    elif line_splitted[0][2:] == 'border_correction':
                        border_correction = eval(line_splitted[1])
                    elif line_splitted[0][2:] == 'FORMAT':
                        fmt_id = line.strip().split(',')[0][-2:]
                        fmt_ds = line.strip().split('"')[1]
                        fmt.append({'ID' : fmt_id, 'DS' : fmt_ds})
                else:
                    col_names = line[1:].strip().split('\t')
            else:
                break
        return CorrelationHeader(ref_mat, tar_mat, sampling_rate, cut_time_before_stim, cut_time_after_stim, binsize, windowsize, border_correction, fmt, col_names)

    def write(self, corr_rec):
        """ Writes one record (one line) in CCG file """
        phases_as_string = ''
        for phase in corr_rec.phases:
            phases_as_string += np.array2string(phase['CH'], max_line_width=1000000) + ':' + str(phase['RS']) + '\t'
        corr_line = '{0}\ttet{1}\t{2}\t{3}\ttet{4}\t{5}\t{6}\t{7}\n'.format(corr_rec.ref_area,
                                                                            corr_rec.ref_tet, 
                                                                            corr_rec.ref_neur+1,
                                                                            corr_rec.tar_area, 
                                                                            corr_rec.tar_tet, 
                                                                            corr_rec.tar_neur+1, 
                                                                            ':'.join([f['ID'] for f in corr_rec.format]), 
                                                                            phases_as_string)
        self.f.write(corr_line)

    def write_header(self):
        self.f.write('##ref_mat={0}\n'.format(self.header.ref_mat))
        self.f.write('##tar_mat={0}\n'.format(self.header.tar_mat))
        self.f.write('##sampling_rate={0}\n'.format(self.header.sampling_rate))
        self.f.write('##cut_time_before_stim={0}\n'.format(self.header.cut_time_before_stim))
        self.f.write('##cut_time_after_stim={0}\n'.format(self.header.cut_time_after_stim))
        self.f.write('##binsize={0}\n'.format(self.header.binsize))
        self.f.write('##windowsize={0}\n'.format(self.header.windowsize))
        self.f.write('##border_correction={0}\n'.format(self.header.border_correction))
        for fmt in self.header.format:
            self.f.write('##FORMAT=<ID={0},Description="{1}">\n'.format(fmt['ID'], fmt['DS']))
        col_names_line = '\t'.join(self.header.col_names)
        self.f.write('#' + col_names_line + '\n')

    def parse_line(self, line):
        """ Parses one record line in CCG file """
        fields = line.split('\t')
        ref_area = fields[0]
        ref_tet = int(fields[1][3:])
        ref_neur = int(fields[2])
        tar_area = fields[3]
        tar_tet = int(fields[4][3:])
        tar_neur = int(fields[5])
        fmt = fields[6].split(':')
        phases_as_string = fields[7:]
        phases = []
        for phase_as_str in phases_as_string:
            phase_dict = dict()
            phase_dict['CH'] = np.fromstring(phase_as_str.split(':')[0][1:-1], sep=' ', dtype=int)
            phase_dict['RS'] = int(phase_as_str.split(':')[1])
            phases.append(phase_dict)
        return CorrelationRecord(ref_area, ref_tet, ref_neur, tar_area, tar_tet, tar_neur, fmt, phases)
        
    def fetch(self):
        """ Get all records from CCG file """
        for line in self.all_lines:
            if not line.startswith('#'):
                yield self.parse_line(line.strip())

            
def list2neo(st1, st2):
    """ Converts two lists of spike times into neo.Spiketrain objects. """
    t_stop = max(max(st1), max(st2)) + 1
    neo_st1 = neo.core.SpikeTrain(st1*pq.us, t_stop=t_stop)
    neo_st2 = neo.core.SpikeTrain(st2*pq.us, t_stop=t_stop)
    return neo_st1, neo_st2


def cch(st1, st2, binsize, windowsize, border_correction):
    """ Cross-correlates two neo spiketrains and returns cross-correlation histogram """
    st2.t_stop = st1.t_stop  # usually the two spiketrains cch() receives
                             # should have the same t_stop but sometimes a
                             # rounding error occured resulting in the two
                             # floating point numbers not being the same
    st1_binned = eph.conversion.BinnedSpikeTrain(st1, binsize=binsize * pq.us)
    st2_binned = eph.conversion.BinnedSpikeTrain(st2, binsize=binsize * pq.us)
    cch_object = eph.spike_train_correlation.cross_correlation_histogram(
        st1_binned, st2_binned, window=[-windowsize, windowsize],
        border_correction=border_correction, binary=False, kernel=None, method='speed')
    cch = [int(val) for val in cch_object[0]]
    return np.array(cch, dtype=int)


def get_cch_for_all_neurons(ref_tet_id, tar_tet_id, ref_spiketimes, tar_spiketimes, phase, binsize, windowsize, border_correction):
    """ Calls cch() function for all combinations of neurouns in the two tetrodes. Each neurons contains data for 4 phases (baseline, 
    study, exp_old, exp_new). Each phase contains 10 spiketrains (1 for each odor (stimulus) presentation). """
    windowsize = windowsize // binsize  # adjusted windowsize
    cch_all_neurons = dict()
    num_ref_spikes_all_neurons = dict()
    for rn in range(len(ref_spiketimes)):
        b = rn+1 if ref_tet_id == tar_tet_id else 0
        for tn in range(b, len(tar_spiketimes)):
            logging.info('CCH @ {0}: RefNeur {1}/{2}, TarNeur {3}/{4}'.format(phase, rn+1, len(ref_spiketimes), tn+1, len(tar_spiketimes)))
            neuron_cch = np.zeros(2*windowsize+1, dtype=int)
            neuron_num_ref_spikes = 0

            for odor in range(len(ref_spiketimes[0])):
                try:
                    min_spike_ts = min(ref_spiketimes[rn][odor][0], tar_spiketimes[tn][odor][0])
                    ref_spiketrain, tar_spiketrain = list2neo(ref_spiketimes[rn][odor] - min_spike_ts, tar_spiketimes[tn][odor] - min_spike_ts)
                    neuron_cch += cch(ref_spiketrain, tar_spiketrain, binsize, windowsize, border_correction)
                    neuron_num_ref_spikes += len(ref_spiketimes[rn][odor])
                except (ValueError, IndexError):  # case when at least one of the two spiketrains contains no neurons at all
                    neuron_cch += np.zeros(2*windowsize+1, dtype=int)
                    neuron_num_ref_spikes += len(ref_spiketimes[rn][odor])
            cch_all_neurons[(rn, tn)] = neuron_cch
            num_ref_spikes_all_neurons[(rn, tn)] = neuron_num_ref_spikes
    return cch_all_neurons, num_ref_spikes_all_neurons


def write_to_ccg(options, P, cch_baseline, cch_study, cch_exp_old, cch_exp_new, ref_spk_baseline, ref_spk_study,
                 ref_spk_exp_old, ref_spk_exp_new):
    """ Writes all cross-correlation histograms to CCG file """
    fmt = [{'ID' : 'CH', 'DS' : 'Cross-Correlation Histogram'}, 
           {'ID' : 'RS', 'DS' : 'Number of spikes in reference spike train'}]
    col_names = ['REFARE', 'REFTET', 'REFNEUR', 'TARARE', 'TARTET', 'TARNEUR', 'FORMAT', 'BASE', 'STUDY', 'EXPOLD', 'EXPNEW']
    ccg_header = CorrelationHeader(options.ref_mat,
                                   options.tar_mat, 
                                   options.sampling_rate, 
                                   options.cut_time_before_stim, 
                                   options.cut_time_after_stim, 
                                   options.binsize,
                                   options.windowsize,
                                   options.border_correction, 
                                   fmt,
                                   col_names)
    ccg_out = CorrelationFile(options.outfile, 'w', header=ccg_header)
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

        corr_rec = CorrelationRecord(P.tetrodes[options.ref_tet_id-1],
                                     options.ref_tet_id, 
                                     neur_pair[0],
                                     P.tetrodes[options.tar_tet_id-1], 
                                     options.tar_tet_id, 
                                     neur_pair[1], 
                                     fmt, 
                                     phases)
        ccg_out.write(corr_rec)

        

        