import matplotlib.pyplot as plt 
import seaborn as sns

def plot_connectivity_heatmap(connectivity_dfs, outfile_root):
    phases = ['Baseline', 'Study', 'Exp_Old_Odors', 'Exp_New_Odors']

    for i in range(len(phases)):
        fig, ax = plt.subplots(figsize=(20, 7))
        sns.heatmap(connectivity_dfs[i], xticklabels=True, cmap='coolwarm', ax=ax, vmin=-1, vmax=1)
        plt.title(phases[i])
        ax.xaxis.tick_top()
        ax.xaxis.set_label_position('top') 
        plt.subplots_adjust(top=0.9)
        plt.ylabel('Neuron [reference tetrode]')
        plt.xlabel('Neuron [target tetrode]')
        plt.savefig(outfile_root + '_' + phases[i] + '.png')
        plt.close