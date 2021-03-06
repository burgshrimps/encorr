import matplotlib.pyplot as plt 
import numpy as np


def plot_cch(corr_rec, ccg_header, conn_rec, ccf_header, workdir):
    plt.rcParams.update(plt.rcParamsDefault)
    plt.rcParams.update({'font.size': 12})
    f, axarr = plt.subplots(1, 4, figsize=(32, 8), sharey=True)
    f.suptitle('Cross-correlation between tetrode ' + str(corr_rec.ref_tet) + ' (' + corr_rec.ref_area + ') and tetrode ' + str(
                corr_rec.tar_tet) + ' (' + corr_rec.tar_area + ')\nNeurons: [' + str(corr_rec.ref_neur) + ', ' + str(
                corr_rec.tar_neur) + ']')

    f.text(0.5, 0.04, 'time lag in ' + str(ccg_header.binsize) + (' us bins'), ha='center')
    f.text(0.08, 0.5, 'Firing probability', va='center', rotation=('vertical'))

    titles = ['Baseline', 'Study Phase', 'Test Phase Old Odors', 'Test Phase New Odors']
    bins = np.array(range(-ccg_header.windowsize, ccg_header.windowsize + 1))
    rects = [[], [], [], []]
    for i in range(len(titles)):
        axarr[i].set_title(titles[i])
        norm_factor = corr_rec.phases[i]['RS'] if corr_rec.phases[i]['RS'] > 0 else 1  # sometimes RS = 0 and then RuntimeWarning occurs, so only divide by number of spikes in ref spiketrain if that number is > 0, otherwise the CCH does not contain any data anyway
        rects[i] = axarr[i].bar(bins, corr_rec.phases[i]['CH'] / norm_factor, color='black', width=1)
        phase_mean = np.mean(corr_rec.phases[i]['CH'] / norm_factor)
        phase_std = np.std(corr_rec.phases[i]['CH'] / norm_factor)

        peak_thr_line = phase_mean + ccf_header.peak_thr * phase_std
        trough_thr_line = phase_mean + ccf_header.trough_thr * phase_std

        axarr[i].plot([-ccg_header.windowsize, ccg_header.windowsize], [peak_thr_line, peak_thr_line],
                       '--', color='limegreen', linewidth=2)
        axarr[i].plot([-ccg_header.windowsize, ccg_header.windowsize], [trough_thr_line, trough_thr_line],
                       '--', color='limegreen', linewidth=2)

        y_lim = axarr[i].get_ylim()[1]
        axarr[i].plot([-ccf_header.center-1, -ccf_header.center-1], [0, y_lim], '--', color='grey')
        axarr[i].plot([ccf_header.center+1, ccf_header.center+1], [0, y_lim], '--', color='grey')

        try:
            for conn in conn_rec.phases[i]:
                if conn['TP'] != '.':
                    if conn['TP'] == 'PK':
                        color = 'orangered'
                    elif conn['TP'] == 'TR':
                        color = 'deepskyblue'
                    bin_idx = conn['BN'] + ccg_header.windowsize
                    axarr[i].text(rects[i][bin_idx].get_x() + rects[i][bin_idx].get_width()/2.0, 1.001 * rects[i][bin_idx].get_height(), '*', color=color, fontsize=25, fontweight='bold', ha='center')
                else:
                    continue
        except AttributeError:
            continue

    _, t = plt.ylim()
    plt.ylim(top=t + (t*0.1))

    plotname = workdir + '/cch_rt' + str(corr_rec.ref_tet) + '_tt' + str(
               corr_rec.tar_tet) + '_rn' + str(corr_rec.ref_neur) + '_tn' + str(corr_rec.tar_neur) + '.png'
    plt.savefig(plotname)
    plt.close()

    


