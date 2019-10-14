import matplotlib.pyplot as plt 
import numpy as np


def plot_cch(corr_rec, ccg_header, conn_rec, ccf_header, workdir):
    plt.rcParams.update(plt.rcParamsDefault)
    plt.rcParams.update({'font.size': 12})
    f, axarr = plt.subplots(1, 4, figsize=(32, 8), sharey=True)
    f.suptitle('Cross-correlation between tetrode ' + str(corr_rec.ref_tet) + ' (' + corr_rec.ref_area + ') and tetrode ' + str(
                corr_rec.tar_tet) + ' (' + corr_rec.tar_area + ')\nNeurons: [' + str(corr_rec.ref_neur) + ', ' + str(
                corr_rec.tar_neur) + ']')

    f.text(0.5, 0.04, 'time lag in ' + str(ccg_header.binsize) + (' ms bins'), ha='center')
    f.text(0.08, 0.5, 'Firing probability', va='center', rotation=('vertical'))

    titles = ['Baseline', 'Study Phase', 'Test Phase Old Odors', 'Test Phase New Odors']
    bins = np.array(range(-ccg_header.windowsize, ccg_header.windowsize + 1))
    rects = [[], [], [], []]
    for i in range(len(titles)):
        axarr[i].set_title(titles[i])
        rects[i] = axarr[i].bar(bins, corr_rec.phases[i]['CH'] / corr_rec.phases[i]['RS'], color='black', width=1)
        
        phase_mean = np.mean(corr_rec.phases[i]['CH'] / corr_rec.phases[i]['RS'])
        phase_std = np.std(corr_rec.phases[i]['CH'] / corr_rec.phases[i]['RS'])
        peak_thr_line = phase_mean + ccf_header.peak_thr * phase_std
        trough_thr_line = phase_mean + ccf_header.trough_thr * phase_std

        axarr[i].plot([-ccg_header.windowsize, ccg_header.windowsize], [peak_thr_line, peak_thr_line],
                       '--', color='limegreen', linewidth=2)
        axarr[i].plot([-ccg_header.windowsize, ccg_header.windowsize], [trough_thr_line, trough_thr_line],
                       '--', color='limegreen', linewidth=2)

        y_lim = axarr[i].get_ylim()[1]
        axarr[i].plot([-6, -6], [0, y_lim], '--', color='grey')
        axarr[i].plot([6, 6], [0, y_lim], '--', color='grey')

    _, t = plt.ylim()
    plt.ylim(top=t + (t*0.1))

    plotname = workdir + '/cch_rt' + str(corr_rec.ref_tet) + '_rn' + str(corr_rec.ref_neur) + '_tt' + str(
               corr_rec.tar_tet) + '_tn' + str(corr_rec.tar_neur) + '.png'
    
    plt.savefig(plotname)

    


