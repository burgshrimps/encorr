import os
import pandas as pd
from scipy.io import savemat
import numpy as np
import csv
import seaborn as sns 
import matplotlib.pyplot as plt

from ENCORR_call import ConnectionFile

plt.style.use('ggplot')

def diag_to_nan(m):
    l = m.shape[0]
    for i in range(l):
        m[i,i,:] = [np.nan, np.nan, np.nan, np.nan]
    return m


def get_intensity_entries(phases):
    intensity_entries = []
    for p in phases:
        if len(p) == 1:
            if p[0]['IN'] == '.':
                intensity_entries.append(0)
            else: 
                intensity_entries.append(p[0]['IN'])
        else:
            peak_present = False
            for conn in p:
                if conn['TP'] == 'PK':
                    peak_present = True
                    break
            if peak_present:
                intensity_entries.append(max([c['IN'] for c in p]))
            else:
                intensity_entries.append(min([c['IN'] for c in p]))
    return intensity_entries

                
def corr_matrix_from_ccf(indir, tet_info):
    neur_count_cum = np.concatenate(([0], np.cumsum(list(tet_info['neur_count']))))  # cumulative neuron counts
    corr_matrix = np.zeros((neur_count_cum[-1], neur_count_cum[-1], 4), dtype=float)
    corr_matrix = diag_to_nan(corr_matrix)
    for file in os.listdir(indir):
        if file.endswith('.ccf'):
            ccf = ConnectionFile(indir + '/' + file, 'r')
            for rec in ccf.fetch():
                i = neur_count_cum[rec.ref_tet-1] + rec.ref_neur - 1
                j = neur_count_cum[rec.tar_tet-1] + rec.tar_neur - 1
                corr_matrix[i,j,:] = get_intensity_entries(rec.phases)
    return corr_matrix, neur_count_cum

def count_conn_per_area(corr_matrix, tet_info, neur_count_cum):
    areas = ['pCA3', 'pCA3x', 'dCA3', 'CA1b', 'pCA1', 'pCA1b', 'dCA1', 'sub', 'cortex', 'x']
    counts = dict()
    npairs = dict()
    for i in range(len(areas)):
        for j in range(i, len(areas)):
            counts[areas[i] + ':' + areas[j]] = np.array([0, 0, 0, 0])
            npairs[areas[i] + ':' + areas[j]] = 0
    
    for rn in range(neur_count_cum[-1]):
        for tn in range(rn+1, neur_count_cum[-1]):
            rt = np.where(rn<neur_count_cum)[0][0]
            tt = np.where(tn<neur_count_cum)[0][0]
            ra = tet_info['area'][rt-1]
            ta = tet_info['area'][tt-1]
            phases_with_conn = np.nonzero(corr_matrix[rn,tn,:])
            try :
                key = ra + ':' + ta
                npairs[key] += 1
                counts[key][phases_with_conn] += 1
            except KeyError:
                key = ta + ':' + ra
                npairs[key] += 1
                counts[key][phases_with_conn] += 1
    return counts, npairs


def to_matlab(tet_info, neur_count_cum, corr_matrix, out_mat, npairs, counts):
    variables = {'neur_corr' : [], 'conn_counts' : []}
    for neur in range(neur_count_cum[-1]):
        tet = np.where(neur<neur_count_cum)[0][0]
        s1 = np.transpose(corr_matrix[:neur,neur,:])
        s2 = np.transpose(corr_matrix[neur,neur:,:])
        variables['neur_corr'].append({'tetrode' : tet,
                                       'neuron' : neur + 1 - neur_count_cum[tet-1],
                                       'area' : tet_info['area'][tet-1],
                                       'xcorr' : np.hstack((s1,s2))})

    for key in npairs:
        norm_npairs = npairs[key] if npairs[key] > 0 else 1
        variables['conn_counts'].append({'pair' : key,
                                         'neur_pairs' : npairs[key],
                                         'counts' : np.vstack((counts[key], counts[key] / norm_npairs))})


    savemat(out_mat, variables)


def node_list_to_csv(tet_info, outfile):
    writer = csv.writer(open(outfile, 'w'))
    header = ['id', 'tetrode', 'neuron', 'area']
    writer.writerow(header)
    i = 1
    for tet in range(len(tet_info['neur_count'])):
        for neuron in range(tet_info['neur_count'][tet]):
            line = [i, tet+1, neuron+1, tet_info['area'][tet]]
            writer.writerow(line)
            i += 1


def edge_list_to_csv(corr_matrix, phase, outfile):
    writer = csv.writer(open(outfile, 'w'))
    header = ['ref', 'tar', 'type', 'weight']
    writer.writerow(header)
    for rn in range(len(corr_matrix)):
        for tn in range(len(corr_matrix[0])):
            curr_entry = corr_matrix[rn, tn, phase]
            if curr_entry > 0:
                writer.writerow([rn+1, tn+1, 'exc', curr_entry])
            elif curr_entry < 0:
                writer.writerow([rn+1, tn+1, 'inh', curr_entry])


def to_csv(corr_matrix, tet_info, out_root):
    node_list_to_csv(tet_info, out_root + '_neurons.csv')
    phases = ['baseline', 'study', 'exp_old', 'exp_new']
    for i in range(4):
        edge_list_to_csv(corr_matrix, i, out_root + '_' + phases[i] + '.csv')


def plot_heatmaps(corr_matrix, out_root):
    phases = ['baseline', 'study', 'exp_old', 'exp_new']
    for i in range(4):
        plt.figure(figsize=(20,8))
        mat = corr_matrix[:,:,i] + corr_matrix[:,:,i].T - np.diag(np.diag(corr_matrix[:,:,i]))
        sns.heatmap(mat, vmin = -1, vmax=1, cmap='coolwarm')
        plt.title(phases[i])
        plt.savefig(out_root + '_heatmap_' + phases[i] + '.png')
        plt.close()


#def plot_pca(corr_matrix, out_root)