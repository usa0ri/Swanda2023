
import pandas as pd
import pysam
import time
import numpy as np
import re
import openpyxl
from openpyxl.styles import numbers,PatternFill
import string

import myRiboSeq.mylib as my

# thresholds for total counts of transcripts
thresholds = [np.inf,64,32,16,8,0]

'''check for NH tags to remove reads that aligned more than once'''
def _check_unique_map(read):
    flag = True
    if read.has_tag('NH'):
        if read.get_tag("NH")>1:
            flag = False
    return flag

'''check for flag to remove reads that aligned as sense (0) or antisense (16)'''
def _check_sense(read):
    return not read.is_reverse

def _calc_mismatch(read,cnt_threshold):
    if read.has_tag('nM'):
        return read.get_tag('nM') < cnt_threshold
    elif read.has_tag('NM'):
        return read.get_tag('NM') < cnt_threshold
    
'''preparation of data for the downstream analyses'''
def prep_data(
    save_dir,
    ref_dir,
    data_dir,
    sp):

    ref = my.Ref(
        data_dir=data_dir,
        ref_dir=ref_dir,
        sp=sp)
    
    if save_dir is None:
        return ref
        
    save_dir = save_dir / 'prep_data'
    if not save_dir.exists():
        save_dir.mkdir()
    
    df_sum = pd.DataFrame(
        np.zeros((len(ref.exp_metadata.df_metadata),2+len(thresholds))),
        columns=['total_reads','filtered_reads','total_tr']+[f'thres{x}_tr' for x in thresholds[1:]],
        index=ref.exp_metadata.df_metadata['sample_name']
        )

    for a in ref.exp_metadata.df_metadata['align_files']:
        smpl_name = ref.exp_metadata.df_metadata.query(f'align_files == "{a}"').sample_name.iloc[-1]
        if (save_dir / f'df_rm_{smpl_name}.csv.gz').exists():
            continue
        print(f'preparing data for {smpl_name}...')
        infile = pysam.AlignmentFile(ref.data_dir / a)

        pointer = infile.tell()
        n_reads_total = infile.count(until_eof=True)
        infile.seek(pointer)

        t0 = time.time()
        flags = [
            read.reference_name in ref.annot.annot_dict.keys() and read.reference_name in ref.id.dict_name.keys()
            for read in infile.fetch(until_eof=True)]
        print(f'Filtering: {np.sum(flags)} left out of {n_reads_total}...')
        
        df_sum.loc[smpl_name,'total_reads'] = n_reads_total
        df_sum.loc[smpl_name,'filtered_reads'] = np.sum(~np.array(flags))

        infile.seek(pointer)
        data_rm = [ [
            read.query_name,
            read.reference_name,
            read.is_reverse,
            read.reference_start,
            read.reference_end,
            read.query_length,
            read.query_alignment_length,
            read.cigarstring,
            read.get_tag('NH'),
            read.get_tag('NM'),
            read.get_tag('nM')]
        for i,read in enumerate(infile.fetch(until_eof=True)) if not flags[i]
        ]
        pd.DataFrame(data_rm,columns=[
            "query_name",
            "tr_id",
            "strand",
            "cut5",
            "cut3",
            "query_length",
            "query_alignment_length",
            "cigarstring",
            "NH","NM","nM"
            ]).to_csv(save_dir / f'df_rm_{smpl_name}.csv.gz',compression="gzip")
        infile.seek(pointer)
        data = [ [
            read.reference_start,# start position of the read (5 end)
            read.reference_end,# end position of the read (3 end)
            read.query_alignment_length,# length of the aligned part of the read
            read.query_length,# read length
            read.is_reverse,# True:+, False:-
            ref.annot.annot_dict[read.reference_name]["start"],# start position of the transcript
            ref.annot.annot_dict[read.reference_name]["stop"]-3# end position of the transcript
            ] for i,read in enumerate(infile.fetch(until_eof=True)) if flags[i]]
        infile.seek(pointer)
            
        attr_names = ["tr_id","seq","read_seq"]
        data_str = [ [
            read.reference_name,
            read.query_alignment_sequence,
            read.query_sequence
            ] for i,read in enumerate(infile.fetch(until_eof=True)) if flags[i]]
        data = np.array(data)
        data = np.hstack((
            data, 
            (data[:, 0] - data[:, 5])[:,np.newaxis],# "dist5_start"
            (data[:, 1] - data[:, 5] -1)[:,np.newaxis],# "dist3_start"
            (data[:, 0] - data[:, 6])[:,np.newaxis],# "dist5_stop"
            (data[:, 1] - data[:, 6] -1)[:,np.newaxis]# "dist3_stop"
        ))
        data = np.hstack((
            data,
            (data[:, 7] % 3)[:,np.newaxis],# frame5
            (data[:, 8] % 3)[:,np.newaxis]# frame3
        ))
        t1 = time.time()
        print(f'\n {np.round(t1-t0,3)} sec')

        '''dataframe of reads
        + index: reads (filtered)
        '''
        df_data = pd.DataFrame(data,columns=[
            "cut5",
            "cut3",
            "length",
            "read_length",
            "strand",
            "start",
            "stop",
            "dist5_start",
            "dist3_start",
            "dist5_stop",
            "dist3_stop",
            "frame5",
            "frame3"
            ])
        df_data.loc[:,"cds_label"] = np.where(df_data['dist5_start']<=-12,"5UTR","CDS")
        df_data.loc[:,"cds_label"] = np.where(df_data['dist5_stop']>-15,"3UTR",df_data.loc[:,"cds_label"])
        df_data.loc[:,attr_names] = data_str

        outfile_name = save_dir / f'df_{smpl_name}.joblib'
        my._mysave(outfile_name,df_data)
            
        '''dictionary of transcripts
        + keys: Ensembl IDs
        + values:   "name": transcript name
                    "tr_info": dataframe ["cds_len","tr_len","counts(within CDS)"]
                    "df": dataframe of reads (only within len_range)
        '''
        t0 = time.time()
        count_tr = df_data.query('cds_label == "CDS"').loc[:,"tr_id"].value_counts()
        df_sum.loc[smpl_name,'total_tr'] = len(count_tr)

        tr_grp = dict(list(df_data.groupby("tr_id")))

        counts_sum = 0
        for i,j in zip(thresholds[:-1],thresholds[1:]):
            is_used = (count_tr <= i) * (count_tr > j)
            count_tr_now = count_tr[is_used]
            print(f'{np.sum(is_used)}/{len(is_used)} transcripts had CDS count > {j} and <= {i}')
            df_sum.loc[smpl_name,f'thres{j}_tr'] = np.sum(is_used)
            counts_sum += count_tr_now.sum()

            dict_tr_now = {
                tr_id:{
                    'name':ref.id.dict_name[tr_id]['symbol'],
                    'tr_info':ref.annot.annot_dict[tr_id],
                    'counts':count_tr_now.iloc[i],
                    'counts_sum':counts_sum,
                    'df':tr_grp[tr_id]}
                for i, tr_id in enumerate(count_tr_now.index)
            }
            outfile_name = save_dir / f'dict_tr_{smpl_name}_tr{j}.joblib'
            my._mysave(outfile_name,dict_tr_now)
            print(f'\n{outfile_name} saved')
        t1 = time.time()
        print(f'\n {np.round(t1-t0,3)} sec')
        print("hoge")
    df_sum.to_csv(save_dir / f'summary.csv.gz',compression='gzip')

    return ref
