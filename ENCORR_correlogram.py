import matplotlib.pyplot as plt 


def plot_cch(corr_rec, header, conn_rec, workdir):
    f, axarr = plt.subplots(1, 4, figsize=(32, 8), sharey=True)
    f.suptitle('Cross-correlation between tetrode ' + str(corr_rec.ref_tet) + ' and tetrode ' + str(
                corr_rec.tar_tet) + '\nNeurons: [' + str(corr_rec.ref_neur) + ', ' + str(corr_rec.tar_neur) + ']')

    f.text(0.5, 0.04, 'time lag in ' + str(BINSIZE) + (' ms bins'), ha='center')
    f.text(0.08, 0.5, 'Firing probability', va='center', rotation=('vertical'))

