# Author:   Nico Alavi, <nico.alavi aT fu-berlin.de>
#           Leibniz-Institute for Neurobiology, Magdeburg
#           Department Functional Architecture of Memory


import logging
import sys
import numpy as np 

from ENCORR_input_parsing import parse_arguments
from ENCORR_load_data import loadparams, loadtet
from ENCORR_cut_out import get_stoi
from ENCORR_cross_correlate import CorrelationFile, get_cch_for_all_neurons, write_to_ccg
from ENCORR_call import ConnectionFile, ConnectionRecord, get_candidates, call_peaks, call_troughs, create_phase_records_firing_prob


def main():
    # Fetch command line argumnets
    options = parse_arguments()
    if not options.sub:
        print('Choose between the modes ("correlate", "call", "correlogram", "connections", "heatmap" or "network").')
        return

    # Set up logging
    logFormatter = logging.Formatter('%(asctime)s %(message)s')
    rootLogger = logging.getLogger()
    rootLogger.setLevel(logging.INFO)

    consoleHandler = logging.StreamHandler()
    consoleHandler.setFormatter(logFormatter)
    rootLogger.addHandler(consoleHandler)

    logging.info('############### Start ENCORR ###############')
    logging.info('CMD: python3 {0}'.format(' '.join(sys.argv)))
    for arg in vars(options):
        logging.info('PARAMETER: {0}, VALUE: {1}'.format(arg, getattr(options, arg)))

    if options.sub == 'correlate':
        logging.info('MODE: correlate')
        logging.info('INDIR: {0}'.format(options.indir))
        logging.info('SAMPLING RATE: {0} kHz'.format(options.sampling_rate))
        logging.info('REFERENCE TETRODE: {0}'.format(options.ref_tet_id))
        logging.info('TARGET TETRODE: {0}'.format(options.tar_tet_id))
        logging.info('OUTPUT FILE: {0}'.format(options.outfile))

        logging.info('# Load parameters')
        P = loadparams(options.indir, options.cut_time_before_stim, options.cut_time_exp_after_resp)
        P.cut_time_before_stim = options.cut_time_before_stim
        P.cut_time_exp_after_resp = options.cut_time_exp_after_resp

        logging.info('# Load tetrodes')
        ref_tet = loadtet(options.indir, options.ref_tet_id, options.sampling_rate)
        tar_tet = loadtet(options.indir, options.tar_tet_id, options.sampling_rate)

        logging.info('# Cut out spiketrains of interest')
        ref_spiketimes_baseline = get_stoi(ref_tet, P, 'baseline')#[:4]
        ref_spiketimes_study = get_stoi(ref_tet, P, 'study')#[:4]
        ref_spiketimes_exp_old = get_stoi(ref_tet, P, 'exp_old')#[:4]
        ref_spiketimes_exp_new = get_stoi(ref_tet, P, 'exp_new')#[:4]
        tar_spiketimes_baseline = get_stoi(tar_tet, P, 'baseline')#[:4]
        tar_spiketimes_study = get_stoi(tar_tet, P, 'study')#[:4]
        tar_spiketimes_exp_old = get_stoi(tar_tet, P, 'exp_old')#[:4]
        tar_spiketimes_exp_new = get_stoi(tar_tet, P, 'exp_new')#[:4]

        logging.info('# Cross-correlate spiketrains of interest')
        cch_baseline, num_ref_spikes_baseline = get_cch_for_all_neurons(options.ref_tet_id, options.tar_tet_id, ref_spiketimes_baseline, tar_spiketimes_baseline, 'baseline',
                                                                        options.binsize, options.windowsize, options.border_correction)
        cch_study, num_ref_spikes_study = get_cch_for_all_neurons(options.ref_tet_id, options.tar_tet_id, ref_spiketimes_study, tar_spiketimes_study, 'study',
                                                                  options.binsize, options.windowsize, options.border_correction)
        cch_exp_old, num_ref_spikes_exp_old = get_cch_for_all_neurons(options.ref_tet_id, options.tar_tet_id, ref_spiketimes_exp_old, tar_spiketimes_exp_old, 'exp_old',
                                                                      options.binsize, options.windowsize, options.border_correction)
        cch_exp_new, num_ref_spikes_exp_new = get_cch_for_all_neurons(options.ref_tet_id, options.tar_tet_id, ref_spiketimes_exp_new, tar_spiketimes_exp_new, 'exp_new',
                                                                      options.binsize, options.windowsize, options.border_correction)

        logging.info('# Write CCH to file')
        write_to_ccg(options.ref_tet_id, options.tar_tet_id, cch_baseline, cch_study, cch_exp_old, cch_exp_new,
                     num_ref_spikes_baseline, num_ref_spikes_study, num_ref_spikes_exp_old, num_ref_spikes_exp_new, options.outfile)

    if options.sub == 'call':
        logging.info('MODE: call')
        logging.info('INPUT CCG FILE: {0}'.format(options.ccg))
        logging.info('OUTPUT CCF FILE: {0}'.format(options.outfile))
        logging.info('PEAK THRESHOLD: {0}'.format(options.peak_thr))
        logging.info('TROUGH THRESHOLD: {0}'.format(options.trough_thr))
        logging.info('PEAK MIN SPIKES: {0}'.format(options.peak_min_spikes))
        logging.info('TROUGH NEIGHBOURS MIN SPIKES: {0}'.format(options.trough_neighbours_min_spikes))
        logging.info('CENTER RANGE: +-{0} BINS'.format(options.center))
        ccg_in = CorrelationFile(options.ccg, 'r')
        ccf_out = ConnectionFile(options.outfile, 'w')
        ccf_out.write_header(options.ccg, options.peak_thr, options.trough_thr, options.peak_min_spikes, options.trough_neighbours_min_spikes, options.center)
        logging.info('# Call Connections')
        for rec in ccg_in.fetch():
            cch_all_phases = np.zeros(len(rec.phases[0]['CH']), dtype=int)
            cch_all_phases_norm = np.zeros(len(rec.phases[0]['CH']), dtype=float)
            for phase in rec.phases:
                cch_all_phases = np.vstack((cch_all_phases, phase['CH']))
                cch_all_phases_norm = np.vstack((cch_all_phases_norm, phase['CH'] / phase['RS']))
            cch_all_phases = cch_all_phases[1:,:]
            cch_all_phases_norm = cch_all_phases_norm[1:,:]

            peak_candidates, trough_candidates = get_candidates(cch_all_phases_norm, options.peak_thr, options.trough_thr)
            peaks_idx = call_peaks(cch_all_phases, peak_candidates, options.peak_min_spikes, options.center)
            troughs_idx = call_troughs(cch_all_phases, trough_candidates, options.trough_neighbours_min_spikes, options.center)

            phases_rec = create_phase_records_firing_prob(cch_all_phases_norm, peaks_idx, troughs_idx)
            conn_rec = ConnectionRecord(rec.ref_tet, rec.ref_neur, rec.tar_tet, rec.tar_neur, 'TP:BN:IN', phases_rec)
            if conn_rec.phases != [[], [], [], []]:  
                ccf_out.write(conn_rec)
        
    #if options.sub == 'correlogram':






main()

