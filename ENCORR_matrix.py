import pandas as pd 
import numpy as np 
from sortedcontainers import SortedDict, SortedSet


def create_sorted_entry_list(entry_set):
    entries_sorted = []
    idx_dict = {}
    idx = 0
    for t in entry_set:
        for n in entry_set[t]:
            entry_id = 'tet' + str(t) + ':' + str(n)
            entries_sorted.append(entry_id)
            idx_dict[entry_id] = idx
            idx += 1

    return entries_sorted, idx_dict



def build_matrices(ccf_in):
    ref_entry_sets = SortedDict()
    tar_entry_sets = SortedDict()
    for rec in ccf_in.fetch():
        try:
            ref_entry_sets[rec.ref_tet].add(rec.ref_neur)
        except KeyError:
            ref_entry_sets[rec.ref_tet] = SortedSet([rec.ref_neur])
        
        try:
            tar_entry_sets[rec.tar_tet].add(rec.tar_neur)
        except KeyError:
            tar_entry_sets[rec.tar_tet] = SortedSet([rec.tar_neur])
    
    ref_entries_sorted, ref_idx = create_sorted_entry_list(ref_entry_sets)
    tar_entries_sorted, tar_idx = create_sorted_entry_list(tar_entry_sets)

    matrices = [np.zeros((len(ref_entries_sorted), len(tar_entries_sorted)), dtype=float) for i in range(4)]
    
    for rec in ccf_in.fetch():
        for i in range(len(rec.phases)):
            for j in range(len(rec.phases[i])):
                if rec.phases[i][j]['TP'] != '.':
                    ref_id = 'tet' + str(rec.ref_tet) + ':' + str(rec.ref_neur)
                    tar_id = 'tet' + str(rec.tar_tet) + ':' + str(rec.tar_neur)
                    matrices[i][ref_idx[ref_id], tar_idx[tar_id]] += rec.phases[i][j]['IN']

    dataframes = []
    for matrix in matrices:
        dataframes.append(pd.DataFrame(matrix, index=ref_entries_sorted, columns=tar_entries_sorted))
    
    return dataframes
    
        
