import argparse
import sys


def parse_arguments(arguments = sys.argv[1:]):
      parser = argparse.ArgumentParser(description='ENCORR')

      subparsers = parser.add_subparsers(help='modes', dest='sub')

      parser_correlate = subparsers.add_parser('correlate',
                                                help='Cross-correlate all neurons in the reference tetrode and all neurons in the target \
                                                      tetrode and saves cross-correlograms in .ccg file.')
      parser_correlate.add_argument('ref_mat',
                                     type=str,
                                     help='MAT file with spiketimes of reference tetrode.')
      parser_correlate.add_argument('tar_mat',
                                     type=str,
                                     help='MAT file with spiketimes of target tetrode.')
      parser_correlate.add_argument('params_mat',
                                     type=str,
                                     help='MAT file with experimental parameters.')
      parser_correlate.add_argument('ref_tet_id',
                                     type=int,
                                     help='ID of reference tetrode.')
      parser_correlate.add_argument('tar_tet_id',
                                     type=int,
                                     help='ID of target tetrode.')                            
      parser_correlate.add_argument('sampling_rate',
                                     type=float,
                                     help='Sampling rate of electrophysiology recording system.')
      parser_correlate.add_argument('outfile',
                                     type=str,
                                     help='Output .ccg file.')
      parser_correlate.add_argument('--cut_time_before_stim',
                                     metavar='INT',
                                     type=int,
                                     default=200,
                                     help='Time in [ms] before stimulus presentation timestamp = start time of each spiketrain of interest.')
      parser_correlate.add_argument('--cut_time_after_stim',
                                     metavar='INT',
                                     type=int,
                                     default=6000,
                                     help='Time in [ms] after stimulus presentation timestamp = end time of each spiketrain of interest.')
      parser_correlate.add_argument('--baseline_end_time',
                                     metavar='INT',
                                     type=int,
                                     default=300000,
                                     help='Time in [ms] before experiment starts. Used to randomly sample baseline STOIs.')
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

      parser_stat = subparsers.add_parser('corr-stat',
                                           help='Plot statistics based on existing CCF file.')
      parser_stat.add_argument('input_dir',
                                type=str,
                                help='Input directory containg CCF files.')
      parser_stat.add_argument('output_dir',
                                type=str,
                                help='Output directory to save .PNG files in.')

      parser_correlogram = subparsers.add_parser('correlogram',
                                                  help='Create correlogram plots from CCG file with marked connections from CCF file.')
      parser_correlogram.add_argument('ccg',
                                       type=str,
                                       help='Input CCG file.')
      parser_correlogram.add_argument('ccf',
                                       type=str,
                                       help='Input CCF file.')
      parser_correlogram.add_argument('workdir',
                                       type=str,
                                       help='Output and working directory.')

      parser_heatmap = subparsers.add_parser('conn-stat',
                                             help='Construct a .mat file containing a field for each neuron with correlations to all other neurons \
                                                   and connection counts per area.')
      parser_heatmap.add_argument('input_dir',
                                  type=str,
                                  help='Directory containing CCF files.')
      parser_heatmap.add_argument('tet_info',
                                  type=str,
                                  help='MAT file containing tetrode information.')
      parser_heatmap.add_argument('out_mat',
                                  type=str,
                                  help='Output MAT file.')

      return parser.parse_args(arguments)
