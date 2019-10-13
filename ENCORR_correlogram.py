import matplotlib.pyplot as plt 


def plot_cch(corr_rec, ccg_header, conn_rec, ccf_header, workdir):
    f, axarr = plt.subplots(1, 4, figsize=(32, 8), sharey=True)
    f.suptitle('Cross-correlation between tetrode ' + str(corr_rec.ref_tet) + ' and tetrode ' + str(
                corr_rec.tar_tet) + '\nNeurons: [' + str(corr_rec.ref_neur) + ', ' + str(corr_rec.tar_neur) + ']')

    f.text(0.5, 0.04, 'time lag in ' + str(ccg_header.binsize) + (' ms bins'), ha='center')
    f.text(0.08, 0.5, 'Firing probability', va='center', rotation=('vertical'))

    titles = ['Baseline', 'Study Phase', 'Test Phase Old Odors', 'Test Phase New Odors']
    rects = [[], [], [], []]
    for i in range(len(titles)):
        axarr[i].set_title(titles[i])
        rects[i] = axarr[i].bar(corr_rec.phases[i]['CH'] / corr_rec.phases[i]['RS'])

    


