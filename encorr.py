# Author:   Nico Alavi, <nico.alavi aT fu-berlin.de>
#           Leibniz-Institute for Neurobiology, Magdeburg
#           Department Functional Architecture of Memory


import logging
import sys
import os
import glob
import numpy as np 
import pickle
import pandas as pd
import matplotlib.pyplot as plt 
import seaborn as sns

from ENCORR_input_parsing import parse_arguments
from ENCORR_load_data import loadparams, loadtet, load_tet_info
from ENCORR_cut_out import get_stoi, get_stoi_baseline_spkcount
from ENCORR_cross_correlate import CorrelationFile, get_cch_for_all_neurons, write_to_ccg
from ENCORR_call import ConnectionFile, ConnectionHeader, ConnectionRecord, get_candidates, call_peaks, call_troughs, create_phase_records
from ENCORR_corr_stat import plot_stat_intensity, plot_stat_bin
from ENCORR_correlogram import plot_cch
from ENCORR_conn_stat import corr_matrix_from_ccf, to_matlab, count_conn_per_area, to_csv, plot_heatmaps, plot_pca, create_summary_csv



def main():
    # Fetch command line argumnets
    options = parse_arguments()
    if not options.sub:
        print('Choose between the modes ("correlate", "call", "correlogram", "corr-stat", "conn-stat").')
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

    if options.sub == 'explore':
        P = loadparams(options.params_mat, options.cut_time_before_stim, options.cut_time_after_stim, 
                       options.baseline_end_time)
        tet_list = glob.glob(options.tet_dir + '/*tet[0-9]*.mat')
        tet_list.sort()
        df = pd.DataFrame()
        i = 0
        axs_idx = [(i,j) for i in range(4) for j in range(4)]
        fig, axs = plt.subplots(4,4, sharex=True, figsize=(20,10))
        for tet_fname in tet_list:
            tet_id = tet_fname[-6:-4] if tet_fname[-6] != 't' else tet_fname[-5:-4] # all samples except LE46
            #tet_id = tet_fname[-12] use for LE46
            tet = loadtet(tet_fname, options.sampling_rate)
            _, num_spikes_baseline = get_stoi(tet, P, 'baseline')
            _, num_spikes_study = get_stoi(tet, P, 'study')
            _, num_spikes_exp_old = get_stoi(tet, P, 'exp_old')
            _, num_spikes_exp_new = get_stoi(tet, P, 'exp_new')

            df_baseline = pd.DataFrame(num_spikes_baseline)
            df_baseline['phase'] = 'baseline'
            df_baseline['tet'] = i
            df_study = pd.DataFrame(num_spikes_study)
            df_study['phase'] = 'study'
            df_study['tet'] = i
            df_exp_old = pd.DataFrame(num_spikes_exp_old)
            df_exp_old['phase'] = 'exp_old'
            df_exp_old['tet'] = i
            df_exp_new= pd.DataFrame(num_spikes_exp_new)
            df_exp_new['phase'] = 'exp_new'
            df_exp_new['tet'] = i

            df = df_baseline.append([df_study, df_exp_old, df_exp_new], ignore_index=True)
            
            df = df.rename(columns={0 : 'num_spikes'})
            sns.boxplot(x='phase', y='num_spikes', data=df, showfliers = False, ax=axs[axs_idx[i]])
            axs[axs_idx[i]].set_ylabel('')    
            axs[axs_idx[i]].set_xlabel('')
            axs[axs_idx[i]].set_title('Tetrode ' + str(tet_id))
            i += 1
        fig.text(0.5, 0.04, 'Phase', ha='center', va='center', fontsize=18)
        fig.text(0.06, 0.5, 'Number of  Spikes in STOIs', ha='center', va='center', rotation='vertical', fontsize=18)
        plt.suptitle(options.name, fontsize=22)
        plt.savefig(options.outfile)

    if options.sub == 'correlate':
        logging.info('MODE: correlate')
        logging.info('REFERENCE TETRODE: {0}'.format(options.ref_mat))
        logging.info('TARGET TETRODE: {0}'.format(options.tar_mat))
        logging.info('PARAMETERS: {0}'.format(options.params_mat))
        logging.info('REF TET ID: {0}'.format(options.ref_tet_id))
        logging.info('TAR TET ID: {0}'.format(options.tar_tet_id))
        logging.info('SAMPLING RATE: {0} kHz'.format(options.sampling_rate))
        logging.info('OUTPUT FILE: {0}'.format(options.outfile))
        logging.info('CUT TIME BEFORE STIM: {0} us'.format(options.cut_time_before_stim))
        logging.info('CUT TIME AFTER STIM: {0} us'.format(options.cut_time_after_stim))
        logging.info('BINSIZE: {0} us'.format(options.binsize))
        logging.info('WINDOWSIZE: {0} us'.format(options.windowsize))
        logging.info('BORDER CORRECTION: {0}'.format(options.border_correction))

        logging.info('# Load parameters')
        P = loadparams(options.params_mat, options.cut_time_before_stim, options.cut_time_after_stim, 
                       options.baseline_end_time)

        logging.info('# Load tetrodes')
        ref_tet = loadtet(options.ref_mat, options.sampling_rate)
        tar_tet = loadtet(options.tar_mat, options.sampling_rate)

        logging.info('# Cut out spiketrains of interest')
        ref_spiketimes_study, ref_spikecount_study = get_stoi(ref_tet, P, 'study')
        ref_spiketimes_exp_old, ref_spikecount_exp_old = get_stoi(ref_tet, P, 'exp_old')
        ref_spiketimes_exp_new, ref_spikecount_exp_new = get_stoi(ref_tet, P, 'exp_new')
        tar_spiketimes_study, tar_spikecount_study = get_stoi(tar_tet, P, 'study')
        tar_spiketimes_exp_old, tar_spikecount_exp_old = get_stoi(tar_tet, P, 'exp_old')
        tar_spiketimes_exp_new, tar_spikecount_exp_new = get_stoi(tar_tet, P, 'exp_new')
        if not options.baseline_spk_count:
            ref_spiketimes_baseline, _ = get_stoi(ref_tet, P, 'baseline')
            tar_spiketimes_baseline, _ = get_stoi(tar_tet, P, 'baseline')
        else:
            ref_spiketimes_baseline = get_stoi_baseline_spkcount(ref_tet, P, ref_spikecount_exp_old, ref_spikecount_exp_new)
            tar_spiketimes_baseline = get_stoi_baseline_spkcount(tar_tet, P, tar_spikecount_exp_old, tar_spikecount_exp_new)
        
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

        logging.info('# Prepare CCF file')
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
        
    if options.sub == 'corr-stat':
        logging.info('MODE: corr-stat')
        logging.info('INPUT DIR: {0}'.format(options.input_dir))
        logging.info('OUTPUT DIR: {0}'.format(options.output_dir))

        stat_bin = [[], [], [], []]
        stat_intensity = [[], [], [], []]
        ref_tet = set()
        tar_tet = set()

        logging.info('# Parse CCF files')
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

        logging.info('# Plot correlation statistics')
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
        if not os.path.exists(options.workdir + '/significant'):
            os.makedirs(options.workdir + '/significant')
        if not os.path.exists(options.workdir + '/not_significant'):
            os.makedirs(options.workdir + '/not_significant')

        for corr_rec in ccg_in.fetch():
            try:
                conn_rec = connections[(corr_rec.ref_tet, corr_rec.ref_neur, corr_rec.tar_tet, corr_rec.tar_neur)]
            except KeyError:
                conn_rec = None
            contains_conn = False
            for phase in conn_rec.phases:
                if phase[0]['TP'] != '.':
                    contains_conn = True
            if contains_conn:  # only plot correlograms with neural connections
                plot_cch(corr_rec, ccg_in.header, conn_rec, ccf_in.header, options.workdir + '/significant')
            else:
                plot_cch(corr_rec, ccg_in.header, conn_rec, ccf_in.header, options.workdir + '/not_significant')
    
    if options.sub == 'conn-stat':
        logging.info('MODE: conn-stat')
        logging.info('CCG INPUT DIR: {0}'.format(options.ccg_input_dir))
        logging.info('CCF INPUT DIR: {0}'.format(options.ccf_input_dir))
        logging.info('TET INFO MAT: {0}'.format(options.tet_info))
        logging.info('OUT ROOT: {0}'.format(options.out_root))
        
        logging.info('Compute connection counts')
        tet_info = load_tet_info(options.tet_info)
        corr_matrix, neur_count_cum = corr_matrix_from_ccf(options.ccf_input_dir, tet_info)
        counts, npairs = count_conn_per_area(corr_matrix, tet_info, neur_count_cum)
        to_matlab(tet_info, neur_count_cum, corr_matrix, options.out_root + '_conn_stat.mat', npairs, counts)

        logging.info('Plot heatmap and PCA')
        id_to_area = to_csv(corr_matrix, tet_info, options.out_root)
        plot_heatmaps(corr_matrix, options.out_root)
        plot_pca(corr_matrix, id_to_area, options.out_root)
        
        create_summary_csv(options.ccf_input_dir, options.ccg_input_dir, options.out_root + '_ccg_wings_stat.tsv')


main()

