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

      parser_call = subparsers.add_parser('call',
                                           help='Call features indicating connectivity between two neurons such as peaks, troughs, \
                                                 and wide troughs from a given cross-correlation histogram.')
      parser_call.add_argument('ccg',
                                type=str,
                                help='Input .ccg file containing cross-correlation histograms.')
      parser_call.add_argument('outfile',
                                type=str,
                                help='Output .ccf file.')
      parser_call.add_argument('--peak_thr',
                                metavar='INT',
                                type=int,
                                default=4,
                                help='Number of standard deviations above the mean the firing probability has to be to be considered \
                                      for peak calling')
      parser_call.add_argument('--trough_thr',
                                metavar='INT',
                                type=int,
                                default=-3,
                                help='Number of standard deviations below the mean the firing probability has to be to be considered \
                                      for trough calling')
      parser_call.add_argument('--peak_min_spikes',
                                metavar='INT',
                                type=int,
                                default=3,
                                help='Number of spikes a peak has to contain at least.')
      parser_call.add_argument('--trough_neighbours_min_spikes',
                                metavar='INT',
                                type=int,
                                default=3,
                                help='Number of spikes a the +- 2 neighbouring bins of a trough have to contain at least.')
      parser_call.add_argument('--center',
                                metavar='INT',
                                type=int,
                                default=5,
                                help='Number of bins around the center bin of the CCH to be considered the region of interest.')
      

      return parser.parse_args(arguments)
