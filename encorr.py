# Author:   Nico Alavi, <nico.alavi aT fu-berlin.de>
#           Leibniz-Institute for Neurobiology, Magdeburg
#           Department Functional Architecture of Memory


import logging
import sys
import os
import numpy as np 
import pickle
import pandas as pd

from ENCORR_input_parsing import parse_arguments
from ENCORR_load_data import loadparams, loadtet, load_tet_info
from ENCORR_cut_out import get_stoi
from ENCORR_cross_correlate import CorrelationFile, get_cch_for_all_neurons, write_to_ccg
from ENCORR_call import ConnectionFile, ConnectionHeader, ConnectionRecord, get_candidates, call_peaks, call_troughs, create_phase_records
from ENCORR_stat import plot_stat_intensity, plot_stat_bin
from ENCORR_correlogram import plot_cch
from ENCORR_conn_stat import corr_matrix_from_ccf, to_matlab, count_conn_per_area



def main():
    # Fetch command line argumnets
    options = parse_arguments()
    if not options.sub:
        print('Choose between the modes ("correlate", "call", "correlogram", "stat", "matrix", "heatmap" or "network").')
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
        logging.info('REFERENCE TETRODE: {0}'.format(options.ref_mat))
        logging.info('TARGET TETRODE: {0}'.format(options.tar_mat))
        logging.info('PARAMETERS: {0}'.format(options.params_mat))
        logging.info('REF TET ID: {0}'.format(options.ref_tet_id))
        logging.info('TAR TET ID: {0}'.format(options.tar_tet_id))
        logging.info('SAMPLING RATE: {0} kHz'.format(options.sampling_rate))
        logging.info('OUTPUT FILE: {0}'.format(options.outfile))
        logging.info('CUT TIME BEFORE STIM: {0} ms'.format(options.cut_time_before_stim))
        logging.info('CUT TIME EXP AFTER RESP: {0} ms'.format(options.cut_time_exp_after_resp))
        logging.info('BINSIZE: {0} ms'.format(options.binsize))
        logging.info('WINDOWSIZE: {0} ms'.format(options.windowsize))
        logging.info('BORDER CORRECTION: {0}'.format(options.border_correction))

        logging.info('# Load parameters')
        P = loadparams(options.params_mat, options.cut_time_before_stim, options.cut_time_exp_after_resp)
        P.cut_time_before_stim = options.cut_time_before_stim
        P.cut_time_exp_after_resp = options.cut_time_exp_after_resp

        logging.info('# Load tetrodes')
        ref_tet = loadtet(options.ref_mat, options.sampling_rate)
        tar_tet = loadtet(options.tar_mat, options.sampling_rate)

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
        write_to_ccg(options, P, cch_baseline, cch_study, cch_exp_old, cch_exp_new,
                     num_ref_spikes_baseline, num_ref_spikes_study, num_ref_spikes_exp_old, num_ref_spikes_exp_new)

    if options.sub == 'call':
        logging.info('MODE: call')
        logging.info('INPUT CCG FILE: {0}'.format(options.ccg))
        logging.info('OUTPUT CCF FILE: {0}'.format(options.outfile))
        logging.info('PEAK THRESHOLD: {0}'.format(options.peak_thr))
        logging.info('TROUGH THRESHOLD: {0}'.format(options.trough_thr))
        logging.info('PEAK MIN SPIKES: {0}'.format(options.peak_min_spikes))
        logging.info('TROUGH NEIGHBOURS MIN SPIKES: {0}'.format(options.trough_neighbours_min_spikes))
        logging.info('CENTER RANGE: +-{0} BINS'.format(options.center))

        fmt = [{'ID' : 'TP', 'DS' : 'Connection Type'}, 
               {'ID' : 'BN', 'DS' : 'Correlogram Bin'}, 
               {'ID' : 'IN', 'DS' : 'Intensity'}]
        col_names = ['REFARE', 'REFTET', 'REFNEUR', 'TARARE', 'TARTET', 'TARNEUR', 'FORMAT', 'BASE', 'STUDY', 'EXPOLD', 'EXPNEW']
        ccf_header = ConnectionHeader(options.ccg, 
                                      options.peak_thr, 
                                      options.trough_thr, 
                                      options.peak_min_spikes, 
                                      options.trough_neighbours_min_spikes, 
                                      options.center, 
                                      fmt,
                                      col_names)

        ccg_in = CorrelationFile(options.ccg, 'r')
        ccf_out = ConnectionFile(options.outfile, 'w', header=ccf_header)


        logging.info('# Call Connections')
        for rec in ccg_in.fetch():
            cch_all_phases = np.zeros(len(rec.phases[0]['CH']), dtype=int)
            cch_all_phases_norm = np.zeros(len(rec.phases[0]['CH']), dtype=float)
            for phase in rec.phases:
                cch_all_phases = np.vstack((cch_all_phases, phase['CH']))
                if phase['RS'] > 0:
                    cch_all_phases_norm = np.vstack((cch_all_phases_norm, phase['CH'] / phase['RS']))
                else:  # case no spikes in reference spiketrain, if no spikes then whole CCH 0 anyways
                    cch_all_phases_norm = np.vstack((cch_all_phases_norm, np.zeros(len(rec.phases[0]['CH']), dtype=float)))
            cch_all_phases = cch_all_phases[1:,:]
            cch_all_phases_norm = cch_all_phases_norm[1:,:]

            peak_candidates, trough_candidates = get_candidates(cch_all_phases_norm, options.peak_thr, options.trough_thr)
            peaks_idx = call_peaks(cch_all_phases, peak_candidates, options.peak_min_spikes, options.center)
            troughs_idx = call_troughs(cch_all_phases, trough_candidates, options.trough_neighbours_min_spikes, options.center)

            phases_rec = create_phase_records(cch_all_phases_norm, peaks_idx, troughs_idx)
            conn_rec = ConnectionRecord(rec.ref_area, rec.ref_tet, rec.ref_neur, rec.tar_area, rec.tar_tet, rec.tar_neur, fmt, phases_rec)
            ccf_out.write(conn_rec)
        
    if options.sub == 'stat':
        logging.info('MODE: stat')
        logging.info('INPUT DIR: {0}'.format(options.input_dir))
        logging.info('OUTPUT DIR: {0}'.format(options.output_dir))

        stat_bin = [[], [], [], []]
        stat_intensity = [[], [], [], []]
        ref_tet = set()
        tar_tet = set()

        for file in os.listdir(options.input_dir):
            if file.endswith('.ccf'):
                ccf_in = ConnectionFile(options.input_dir + '/' + file, 'r')
                for rec in ccf_in.fetch():
                    for i in range(len(rec.phases)):
                        for conn in rec.phases[i]:
                            if conn['TP'] != '.':
                                ref_tet.add(rec.ref_tet)
                                tar_tet.add(rec.tar_tet)
                                stat_bin[i].append(conn['BN'])
                                stat_intensity[i].append(conn['IN'])

        plot_stat_intensity(stat_intensity, ref_tet, tar_tet, options.output_dir)
        plot_stat_bin(stat_bin, ref_tet, tar_tet, options.output_dir)

    if options.sub == 'correlogram':
        logging.info('MODE: correlogram')
        logging.info('CCG FILE: {0}'.format(options.ccg))
        logging.info('CCF FILE: {0}'.format(options.ccf))
        logging.info('WORKDIR: {0}'.format(options.workdir))

        ccg_in = CorrelationFile(options.ccg, 'r')
        ccf_in = ConnectionFile(options.ccf, 'r')

        logging.info('# Fetch CCF records')
        connections = dict()
        for conn_rec in ccf_in.fetch():
            connections[(conn_rec.ref_tet, conn_rec.ref_neur, conn_rec.tar_tet, conn_rec.tar_neur)] = conn_rec 

        logging.info('# Plot correlograms')
        if not os.path.exists(options.workdir):
            os.makedirs(options.workdir)

        for corr_rec in ccg_in.fetch():
            try:
                conn_rec = connections[(corr_rec.ref_tet, corr_rec.ref_neur, corr_rec.tar_tet, corr_rec.tar_neur)]
            except KeyError:
                conn_rec = None
            plot_cch(corr_rec, ccg_in.header, conn_rec, ccf_in.header, options.workdir)
    
    if options.sub == 'conn-stat':
        tet_info = load_tet_info(options.tet_info)
        corr_matrix, neur_count_cum = corr_matrix_from_ccf(options.input_dir, tet_info)
        counts, npairs = count_conn_per_area(corr_matrix, tet_info, neur_count_cum)
        to_matlab(tet_info, neur_count_cum, corr_matrix, options.out_mat, npairs, counts)


main()

