import os
import pandas as pd
from scipy.io import savemat
import numpy as np

from ENCORR_call import ConnectionFile


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
    counts = {'pCA3:pCA3' : np.array([0, 0, 0, 0]),
              'pCA3:dCA3' : np.array([0, 0, 0, 0]),
              'pCA3:pCA1' : np.array([0, 0, 0, 0]),
              'pCA3:dCA1' : np.array([0, 0, 0, 0]),
              'dCA3:dCA3' : np.array([0, 0, 0, 0]),
              'dCA3:pCA1' : np.array([0, 0, 0, 0]),
              'dCA3:dCA1' : np.array([0, 0, 0, 0]),
              'pCA1:pCA1' : np.array([0, 0, 0, 0]),
              'pCA1:dCA1' : np.array([0, 0, 0, 0]),
              'dCA1:dCA1' : np.array([0, 0, 0, 0])}

    npairs = {'pCA3:pCA3' : 0,
              'pCA3:dCA3' : 0,
              'pCA3:pCA1' : 0,
              'pCA3:dCA1' : 0,
              'dCA3:dCA3' : 0,
              'dCA3:pCA1' : 0,
              'dCA3:dCA1' : 0,
              'pCA1:pCA1' : 0,
              'pCA1:dCA1' : 0,
              'dCA1:dCA1' : 0}

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