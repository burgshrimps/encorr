import argparse
import sys


def parse_arguments(arguments = sys.argv[1:]):
    parser = argparse.ArgumentParser(description='ENCORR')

    subparsers = parser.add_subparsers(help='modes', dest='sub')

    parser_correlate = subparsers.add_parser('correlate',
                                              help='Cross-correlate all neurons in the reference tetrode and all neurons in the target \
                                                    tetrode and saves cross-correlograms in .ccg file.')
    parser_correlate.add_argument('indir',
                                   type=str,
                                   help='Input directory containing .mat files for experimental parameters and spiketimes.')
    parser_correlate.add_argument('sampling_rate',
                                   type=int,
                                   help='Sampling rate of electrophysiology recording system.')
    parser_correlate.add_argument('ref_tet_id',
                                   type=int,
                                   help='ID of reference tetrode.')
    parser_correlate.add_argument('tar_tet_id',
                                   type=int,
                                   help='ID of target tetrode.')
    parser_correlate.add_argument('outfile',
                                   type=str,
                                   help='Output .ccg file.')
    parser_correlate.add_argument('--cut_time_before_stim',
                                   metavar='INT',
                                   type=int,
                                   default=200,
                                   help='Time in [ms] before stimulus presentation timestamp = start time of each spiketrain of interest.')
    parser_correlate.add_argument('--cut_time_exp_after_resp',
                                   metavar='INT',
                                   type=int,
                                   default=1000,
                                   help='Time in [ms] after response timestamp in exp phase = end time of each spiketrain of interest \
                                         in exp phase.')
    parser_correlate.add_argument('--binsize',
                                   metavar='INT',
                                   type=int,
                                   default=1,
                                   help='Size of bins in [ms] for binning the spiketimes during cross-correlation analysis.')
    parser_correlate.add_argument('--windowsize',
                                   metavar='INT',
                                   type=int,
                                   default=50,
                                   help='Size of shifting window in [ms] for cross-correlation analysis.')
    parser_correlate.add_argument('--border_correction',
                                   metavar='BOOL',
                                   type=bool,
                                   default=False,
                                   help='Whether to correct for the border effect during cross-correlation analysis. See elephant software \
                                         library for more information.')

    return parser.parse_args(arguments)
