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


class ConnectionFile:

    def __init__(self, outfile, read_write):
        self.outfile = outfile
        self.f = open(self.outfile, read_write)

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

    def write_header(self, ccg, peak_thr, trough_thr, peak_min_spikes, trough_neighbours_min_spikes, center):
        self.f.write('##ccg={0}\n'.format(ccg))
        self.f.write('##peak_thr={0}\n'.format(peak_thr))
        self.f.write('##trough_thr={0}\n'.format(trough_thr))
        self.f.write('##peak_min_spikes={0}\n'.format(peak_min_spikes))
        self.f.write('##trough_neighbours_min_spikes={0}\n'.format(trough_neighbours_min_spikes))
        self.f.write('##center={0}\n'.format(center))
        self.f.write('##FORMAT=<ID=TP,Description="Connection Type">\n')
        self.f.write('##FORMAT=<ID=BN,Description="Correlogram Bin">\n')
        self.f.write('##FORMAT=<ID=IN,Description="Intensity">\n')
        self.f.write('#REFTET\tREFNEUR\tTARTET\tTARNEUR\tFORMAT\tBASE\tSTUDY\tEXPOLD\tEXPNEW\n')

    def write(self, conn_rec):
        phases_as_str = self.phases_rec_to_str(conn_rec.phases)
        conn_line = 'tet{0}\t{1}\ttet{2}\t{3}\t{4}\t{5}\n'.format(conn_rec.ref_tet, 
                                                      conn_rec.ref_neur,
                                                      conn_rec.tar_tet,
                                                      conn_rec.tar_neur,
                                                      conn_rec.format,
                                                      phases_as_str)
        self.f.write(conn_line)
    
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


def create_phase_records_firing_prob(cch_all_phases_norm, peak_indices, trough_indices):
    phases = [[], [], [], []]
    window_length = int(len(cch_all_phases_norm[0]) / 2)

    for pidx in peak_indices:
        chance_firing_prob = np.mean(np.concatenate((cch_all_phases_norm[pidx[0],:29], cch_all_phases_norm[pidx[0], 81:])))
        phases[pidx[0]].append({'TP' : 'PK', 
                                'BN' : pidx[1] - window_length,
                                'IN' : np.round(cch_all_phases_norm[pidx] - chance_firing_prob, 3)})
    
    for tidx in trough_indices:
        chance_firing_prob = np.mean(np.concatenate((cch_all_phases_norm[tidx[0],:29], cch_all_phases_norm[tidx[0], 81:])))
        phases[tidx[0]].append({'TP' : 'TR', 
                                'BN' : tidx[1] - window_length,
                                'IN' : np.round(cch_all_phases_norm[tidx] - chance_firing_prob, 3)})

    return phases


        


    


