from pathlib import Path

cur_dir = Path(__file__)
save_dir = cur_dir.parent / "result/res20230521/Swanda2023"
if not save_dir.exists():
    save_dir.mkdir(parents=True)

from myRiboSeq import Swanda2023

Swanda2023.plot_go(
    save_dir=save_dir,
    fpaths = [
        save_dir / 'goatools' / 'increase.xlsx',
        save_dir / 'goatools' / 'decrease.xlsx'
    ],
    num_plot = 10
)


from myRiboSeq import myprep

ref = myprep.prep_data(
    save_dir=None,
    ref_dir = 'ref/Mus_musculus_106_saori',
    data_dir="data/MiniSeq_Qianlab/data_Swanda2023")

smpls = ref.exp_metadata.df_metadata['sample_name'].values

Swanda2023.cpm(
    save_dir=save_dir,
    load_dir=save_dir / 'prep_data',
    ref=ref
)

print("hoge")

pairs = [
    ('CTNS_KD','NT'),
    ('Cysteamine','NT'),
]

Swanda2023.volcano_plot(
    save_dir=save_dir,
    load_dir=save_dir / 'cpm',
    pairs=pairs,
    threshold_fc=2,
    threshold_pval=0.05
)

pairs = [
    ('mRNA_-AA','mRNA_NT'),
    ('mRNA_-Cys','mRNA_NT'),
]

Swanda2023.volcano_plot(
    save_dir=save_dir,
    load_dir=save_dir / 'cpm',
    pairs=pairs,
    threshold_fc=2,
    threshold_pval=0.05
)

print("hoge")


