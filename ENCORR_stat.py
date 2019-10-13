import seaborn as sns
import matplotlib.pyplot as plt 
import numpy as np

def plot_stat_intensity(stat_intensity, ref_tet, tar_tet, workdir):
    
    titles = ['Baseline', 'Study', 'Exp Old', 'Exp New']
    fig = plt.figure(figsize=(20,8))
    fig.subplots_adjust(hspace=0.4, wspace=0.4)
    fig.suptitle('Connection Intensity Distributions\nRefTet={0}, TarTet={1}'.format(ref_tet, tar_tet))
    for i in range(4):
        plt.subplot(221 + i)
        sns.distplot(stat_intensity[i], bins=np.arange(-1, 1.1, 0.1), kde=False, hist_kws={'edgecolor' : 'white'})
        sns.despine()
        plt.title(titles[i] + ' (n={0})'.format(len(stat_intensity[i])))
        plt.xlabel('Connection Intensity')
    plt.savefig(workdir + '/connection_intensity_distribution.png')


def plot_stat_bin(stat_bin, ref_tet, tar_tet, workdir):
    titles = ['Baseline', 'Study', 'Exp Old', 'Exp New']
    fig = plt.figure(figsize=(20,8))
    fig.subplots_adjust(hspace=0.4, wspace=0.4)
    fig.suptitle('Connection Timing Distributions\nRefTet={0}, TarTet={1}'.format(ref_tet, tar_tet))
    for i in range(4):
        plt.subplot(221 + i)
        sns.distplot(stat_bin[i], bins=np.arange(-5.5, 6.5, 1), kde=False, hist_kws={'edgecolor' : 'white'})
        sns.despine()
        plt.title(titles[i] + ' (n={0})'.format(len(stat_bin[i])))
        plt.xlabel('Bin')
    plt.savefig(workdir + '/connection_timing_distribution.png')