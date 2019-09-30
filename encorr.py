# Author:   Nico Alavi, <nico.alavi aT fu-berlin.de>
#           Leibniz-Institute for Neurobiology, Magdeburg
#           Department Functional Architecture of Memory


import logging
import sys
import numpy as np 

from ENCORR_input_parsing import parse_arguments
from ENCORR_load_data import loadparams, loadtet
from ENCORR_cut_out import get_stoi
from ENCORR_cross_correlate import get_cch_for_all_neurons, create_cch_struct


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
        ref_spiketimes_baseline = get_stoi(ref_tet, P, 'baseline')[:4]
        ref_spiketimes_study = get_stoi(ref_tet, P, 'study')[:4]
        ref_spiketimes_exp_old = get_stoi(ref_tet, P, 'exp_old')[:4]
        ref_spiketimes_exp_new = get_stoi(ref_tet, P, 'exp_new')[:4]
        tar_spiketimes_baseline = get_stoi(tar_tet, P, 'baseline')[:4]
        tar_spiketimes_study = get_stoi(tar_tet, P, 'study')[:4]
        tar_spiketimes_exp_old = get_stoi(tar_tet, P, 'exp_old')[:4]
        tar_spiketimes_exp_new = get_stoi(tar_tet, P, 'exp_new')[:4]

        logging.info('# Cross-correlate spiketrains of interest')
        cch_baseline, num_ref_spikes_baseline = get_cch_for_all_neurons(options.ref_tet_id, options.tar_tet_id, ref_spiketimes_baseline, tar_spiketimes_baseline, 'baseline',
                                                                        options.binsize, options.windowsize, options.border_correction)
        cch_study, num_ref_spikes_study = get_cch_for_all_neurons(options.ref_tet_id, options.tar_tet_id, ref_spiketimes_study, tar_spiketimes_study, 'study',
                                                                  options.binsize, options.windowsize, options.border_correction)
        cch_exp_old, num_ref_spikes_exp_old = get_cch_for_all_neurons(options.ref_tet_id, options.tar_tet_id, ref_spiketimes_exp_old, tar_spiketimes_exp_old, 'exp_old',
                                                                      options.binsize, options.windowsize, options.border_correction)
        cch_exp_new, num_ref_spikes_exp_new = get_cch_for_all_neurons(options.ref_tet_id, options.tar_tet_id, ref_spiketimes_exp_new, tar_spiketimes_exp_new, 'exp_new',
                                                                      options.binsize, options.windowsize, options.border_correction)

        logging.info('# Create CCH struct')
        cch_collection = create_cch_struct(options.ref_tet_id, options.tar_tet_id, cch_baseline, cch_study, cch_exp_old, cch_exp_new,
                                           num_ref_spikes_baseline, num_ref_spikes_study, num_ref_spikes_exp_old, num_ref_spikes_exp_new)




main()

