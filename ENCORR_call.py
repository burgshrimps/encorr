import numpy as np
from scipy import stats
import warnings

class ConnectionRecord:

    def __init__(self, ref_tet, ref_neur, tar_tet, tar_neur, fmt, phases):
        self.ref_tet = ref_tet
        self.ref_neur = ref_neur
        self.tar_tet = tar_tet
        self.tar_neur = tar_neur
        self.format = fmt
        self.phases = phases


class ConnectionHeader:

    def __init__(self, ccg, peak_thr, trough_thr, peak_min_spikes, trough_neighbours_min_spikes, center, fmt, col_names):
        self.ccg = ccg
        self.peak_thr = peak_thr
        self.trough_thr = trough_thr
        self.peak_min_spikes = peak_min_spikes
        self.trough_neighbours_min_spikes = trough_neighbours_min_spikes
        self.center = center
        self.format = fmt
        self.col_names = col_names


class ConnectionFile:

    def __init__(self, file, read_write, header=None):
        self.header = header
        self.f = open(file, read_write)

        if read_write == 'w' and self.header != None:
            self.write_header(self.header)

        elif read_write == 'r':
            self.all_lines = self.f.readlines()
            self.header = self.read_header()
            self.f.close()
        
    def phases_rec_to_str(self, phases_rec):
        all_phases_as_str = ''
        for phase in phases_rec:
            if len(phase) == 0:
                all_phases_as_str += '.:.:.\t'
            else:
                phase_str = ''
                for conn_entry in phase:
                    phase_str += conn_entry['TP'] + ':' + str(conn_entry['BN']) + ':' + str(conn_entry['IN']) + ':'
                all_phases_as_str += phase_str[:-1] + '\t'
        return all_phases_as_str[:-1]

    def write_header(self, header):
        self.f.write('##ccg={0}\n'.format(header.ccg))
        self.f.write('##peak_thr={0}\n'.format(header.peak_thr))
        self.f.write('##trough_thr={0}\n'.format(header.trough_thr))
        self.f.write('##peak_min_spikes={0}\n'.format(header.peak_min_spikes))
        self.f.write('##trough_neighbours_min_spikes={0}\n'.format(header.trough_neighbours_min_spikes))
        self.f.write('##center={0}\n'.format(header.center))
        for fmt in header.format:
            self.f.write('##FORMAT=<ID={0},Description="{1}">\n'.format(fmt['ID'], fmt['DS']))
        col_names_line = '\t'.join(header.col_names)
        self.f.write('#' + col_names_line + '\n')

    def write(self, conn_rec):
        phases_as_str = self.phases_rec_to_str(conn_rec.phases)
        conn_line = 'tet{0}\t{1}\ttet{2}\t{3}\t{4}\t{5}\n'.format(conn_rec.ref_tet, 
                                                      conn_rec.ref_neur,
                                                      conn_rec.tar_tet,
                                                      conn_rec.tar_neur,
                                                      ':'.join([f['ID'] for f in conn_rec.format]),
                                                      phases_as_str)
        self.f.write(conn_line)
    
    def read_header(self):
        fmt = []
        for line in self.all_lines:
            if line.startswith('#'):
                if line.startswith('##'):
                    line_splitted = line.strip().split('=')
                    if line_splitted[0][2:] == 'ccg':
                        ccg = line_splitted[1]
                    elif line_splitted[0][2:] == 'peak_thr':
                        peak_thr = int(line_splitted[1])
                    elif line_splitted[0][2:] == 'trough_thr':
                        trough_thr = int(line_splitted[1])
                    elif line_splitted[0][2:] == 'peak_min_spikes':
                        peak_min_spikes = int(line_splitted[1])
                    elif line_splitted[0][2:] == 'trough_neighbours_min_spikes':
                        trough_neighbours_min_spikes = int(line_splitted[1])
                    elif line_splitted[0][2:] == 'center':
                        center = int(line_splitted[1])
                    elif line_splitted[0][2:] == 'FORMAT':
                        fmt_id = line.strip().split(',')[0][-2:]
                        fmt_ds = line.strip().split('"')[1]
                        fmt.append({'ID' : fmt_id, 'DS' : fmt_ds})
                else:
                    col_names = line[1:].strip().split('\t')
            else:
                break
        
        return ConnectionHeader(ccg, peak_thr, trough_thr, peak_min_spikes, trough_neighbours_min_spikes, center, fmt, col_names)

    def parse_line(self, line):
        fields = line.split('\t')
        ref_tet = int(fields[0][3:])
        ref_neur = int(fields[1])
        tar_tet = int(fields[2][3:])
        tar_neur = int(fields[3])
        fmt = fields[4].split(':')
        phases_as_string = fields[5:]
        phases = []
        for phase_as_str in phases_as_string:
            conn_in_phase = []
            phase_as_str_splitted = phase_as_str.split(':')
            idx = np.arange(0, len(phase_as_str_splitted)+1, 3)
            num_conn = int(len(phase_as_str_splitted) / 3)
            for i in range(num_conn):
                try:
                    conn_in_phase.append({'TP' : phase_as_str_splitted[idx[i]:idx[i+1]][0], 
                                          'BN' : int(phase_as_str_splitted[idx[i]:idx[i+1]][1]),
                                          'IN' : float(phase_as_str_splitted[idx[i]:idx[i+1]][2])})
                except ValueError:
                    conn_in_phase.append({'TP' : '.', 
                                          'BN' : '.',
                                          'IN' : '.'})
            phases.append(conn_in_phase)
        return ConnectionRecord(ref_tet, ref_neur, tar_tet, tar_neur, fmt, phases)

    def fetch(self):
        for line in self.all_lines:
            if not line.startswith('#'):
                yield self.parse_line(line.strip())
    

def get_candidates(cch_all_phases_norm, peak_thr, trough_thr):
    warnings.filterwarnings('ignore')
    try:    
        z_score_all_phases = stats.zscore(cch_all_phases_norm, axis=1)
        peak_candidates = np.where(z_score_all_phases > peak_thr)
        trough_candidates = np.where(z_score_all_phases < trough_thr)
    except RuntimeWarning:
        peak_candidates = []
        trough_candidates = []
    return list(zip(peak_candidates[0], peak_candidates[1])), list(zip(trough_candidates[0], trough_candidates[1]))


def call_peaks(cch_all_phases, peak_candidates, min_num_spikes, center):
    num_bins = len(cch_all_phases[0])
    center_start = int(num_bins/2) - center
    center_end = int(num_bins/2) + center

    bins_with_peaks = []
    phase_with_peaks_out_center = set()
    for candidate in peak_candidates:
        if cch_all_phases[candidate] >= min_num_spikes:
            bins_with_peaks.append(candidate)
            if candidate[1] < center_start or candidate[1] > center_end:
                phase_with_peaks_out_center.add(candidate[0])

    bins_with_peaks_final = []
    for peak_bin in bins_with_peaks:
        if peak_bin[0] not in phase_with_peaks_out_center:
            bins_with_peaks_final.append(peak_bin)
    
    return bins_with_peaks_final


def call_troughs(cch_all_phases, trough_candidates, min_num_spikes_neighbour, center):
    num_bins = len(cch_all_phases[0])
    center_start = int(num_bins/2) - center
    center_end = int(num_bins/2) + center
    neighbours = [-2, -1, 1, 2]

    bins_with_troughs = []
    phase_with_troughs_out_center = set()
    for candidate in trough_candidates:
        if candidate[1] + neighbours[0] >= 0 and candidate[1] + neighbours[-1] <= num_bins-1:  # check whether troughs goes beyond borders of CCH
            neighbouring_bins_crit = True
            for i in neighbours:
                if cch_all_phases[candidate[0], candidate[1]+i] < min_num_spikes_neighbour:
                    neighbouring_bins_crit = False
            if neighbouring_bins_crit:
                bins_with_troughs.append(candidate)
                if candidate[1] < center_start or candidate[1] > center_end:
                    phase_with_troughs_out_center.add(candidate[0])

    bins_with_troughs_final = []
    for trough_bin in bins_with_troughs:
        if trough_bin[0] not in phase_with_troughs_out_center:
            bins_with_troughs_final.append(trough_bin)

    return bins_with_troughs_final


def create_phase_records(cch_all_phases_norm, peak_indices, trough_indices):
    phases = [[], [], [], []]
    window_length = int(len(cch_all_phases_norm[0]) / 2)

    for pidx in peak_indices:
        chance_firing_prob = np.mean(np.concatenate((cch_all_phases_norm[pidx[0],:29], cch_all_phases_norm[pidx[0], 81:])))
        phases[pidx[0]].append({'TP' : 'PK', 
                                'BN' : pidx[1] - window_length,
                                'IN' : np.round(-(chance_firing_prob - cch_all_phases_norm[pidx]), 3)})
    
    for tidx in trough_indices:
        chance_firing_prob = np.mean(np.concatenate((cch_all_phases_norm[tidx[0],:29], cch_all_phases_norm[tidx[0], 81:])))
        phases[tidx[0]].append({'TP' : 'TR', 
                                'BN' : tidx[1] - window_length,
                                'IN' : np.round(-(chance_firing_prob - cch_all_phases_norm[tidx]), 3)})

    return phases


        


    


