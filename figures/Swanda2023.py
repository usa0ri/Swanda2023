from pathlib import Path
import pandas as pd
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties

import myRiboSeq.mylib as my

# thresholds for total counts of transcripts
thresholds = [np.inf,64,32,16,8,0]
color_frame = my.color_frame
color_region = my.color_region
codon_table = my.codon2aa_table

font_path = '/usr/share/fonts/truetype/msttcorefonts/Arial.ttf'
font_prop = FontProperties(fname=font_path)
plt.rcParams['font.family'] = font_prop.get_name()
plt.rcParams['text.usetex'] = False
plt.rcParams['mathtext.fontset'] = 'custom'
plt.rcParams['mathtext.it'] = 'Arial:italic'
plt.rcParams['mathtext.rm'] = font_prop.get_name()


'''for edgeR inputs (raw count of each gene)'''
def cpm(
    save_dir,
    load_dir,
    ref
):
    save_dir = save_dir / 'cpm'
    if not save_dir.exists():
        save_dir.mkdir()
    
    df_out = []
    for s in ref.exp_metadata.df_metadata['sample_name']:
        print(f'generating count data matrix in {s}...')
        dict_tr,label_tr = my._load_dict_tr(load_dir,s,0)
        # counts = np.array([
        #     v['counts']
        #     for v in dict_tr.values()
        # ])
        counts = np.array([
            len(v['df'])
            for v in dict_tr.values()
        ])
        tr_lists = list(dict_tr.keys())

        if len(df_out) == 0:
            df_out = pd.DataFrame(counts,index=tr_lists,columns=[s])
        else:
            df_out = pd.merge(
                df_out,
                pd.DataFrame(counts,index=tr_lists,columns=[s]),
                how="inner",left_index=True,right_index=True,
                copy=False
            )
    df_out.iloc[:,:6].to_csv(save_dir / f'count_mat_1.csv.gz')
    df_out.iloc[:,6:].to_csv(save_dir / f'count_mat_2.csv.gz')


def volcano_plot(
    save_dir,
    load_dir,
    pairs,
    threshold_fc,
    threshold_pval
):
    save_dir = save_dir / 'volcano_plot'
    if not save_dir.exists():
        save_dir.mkdir()

    outfile_name = save_dir / f'volcano_plot_{pairs[0][0]}_{pairs[1][0]}.pdf'
    pdf = PdfPages(outfile_name)
    fig, axs = plt.subplots(1,2,figsize=(10,5),sharex=True,sharey=True)
    for i,pair in enumerate(pairs):
        df = pd.read_csv( load_dir / f'exact_test_{pair[1]}_{pair[0]}.csv', header=0, index_col=0 )

        idx_high = df.query(f'(logFC >= {np.log2(threshold_fc)}) and (PValue <= {threshold_pval})').index
        idx_low = df.query(f'(logFC <= {np.log2(1/threshold_fc)}) and (PValue <= {threshold_pval})').index
        idx_medium = np.array([x for x in df.index if (x not in idx_high) and (x not in idx_low)])

        for g,c in zip([idx_medium,idx_high,idx_low],['#808080',"#FF0000","#0000FF"]):
            axs[i].scatter(
                df.loc[g,'logFC'],
                df.loc[g,'PValue'].apply(lambda x: -np.log10(x)),
                color=c
            )
        axs[i].set_xlabel('Log fold change',fontsize=15)
        axs[i].set_ylabel('-log10(P value)',fontsize=15)
        axs[i].set_title(f'{pair[0]} vs {pair[1]}',fontsize=15)
    max_xlim = np.max([np.abs(x) for x in axs[i].get_xlim()])
    axs[i].set_xlim(-max_xlim,max_xlim)
    fig.savefig(pdf,format='pdf')
    pdf.close()


def plot_go(
    save_dir,
    fpaths,
    num_plot
):
    save_dir = save_dir / 'plot_go'
    if not save_dir.exists():
        save_dir.mkdir()
    
    df_plots = []
    for fpath in fpaths:
        df = pd.read_excel(fpath, sheet_name=0,header=0)
        df = df.query('NS == "BP"')

        pop_size = df['ratio_in_pop'].apply(lambda x: int(x.split('/')[0])).values
        study_size = df['ratio_in_study'].apply(lambda x: int(x.split('/')[0])).values
        p_fdr_bh = df['p_fdr_bh'].apply(lambda x: -np.log10(x)).values
        fold_enrich = df['ratio_in_study'].apply(lambda x: int(x.split('/')[0]) / int(x.split('/')[1])).values \
            / df['ratio_in_pop'].apply(lambda x: int(x.split('/')[0]) / int(x.split('/')[1])).values
        go_names = (df['name'] + ' (' + df['GO'] + ')').values

        idx = (fold_enrich > 1) * (p_fdr_bh > -np.log10(0.05))
        if np.all(~idx):
            df_plots.append([])
            continue
        elif np.sum(idx)>=num_plot:
            idx_ = np.where(idx)[0]
            idx = idx_[:num_plot]

        df_plot = pd.DataFrame(
            {
                'Fold enrichment':fold_enrich,
                '-log10(FDR)':p_fdr_bh,
                'Number of genes':study_size,
                'Number of genes in GO': pop_size
            },index=go_names
        ).iloc[idx,:]
        df_plot.sort_values(by='Fold enrichment',inplace=True)

        df_plot.to_csv(save_dir / f'{fpath.stem}_plot.csv.gz')
        df_plots.append(df_plot)
    
    outfile_name = save_dir / 'go_plot.pdf'
    pdf = PdfPages(outfile_name)
    fig,axs = plt.subplots(2,1,figsize=(7.5,5))

    for i,(df_plot,fpath) in enumerate(zip(df_plots,fpaths)):
        if len(df_plot) == 0:
            axs[i].set_visible(False)
            continue
        # bar plots
        df_plot.plot.barh(
            y='Fold enrichment',
            color='#808080',
            ax=axs[i],
            legend=False
        )
        axs[i].set_title(fpath.stem)
        axs[i].set_xlabel('Fold enrichment')
    fig.tight_layout()
    fig.savefig(pdf,format='pdf')
    pdf.close()


