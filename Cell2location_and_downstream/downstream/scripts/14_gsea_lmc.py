
import gseapy as gp
import pandas as pd

# Load ranks file
ranks = pd.read_csv(
    '/path/to/data/LMCs_ranks_prostate.tsv',
    sep='\t'
)

# Read as text and replace /t with actual tab
with open('/path/to/data/LMCs_ranks_prostate.tsv', 'r') as f:
    content = f.read().replace('/t', '\t')

# Parse from the fixed string
from io import StringIO
ranks = pd.read_csv(StringIO(content), sep='\t')

print("Columns:", ranks.columns.tolist())
print("Shape:", ranks.shape)
print("\nHead:")
print(ranks.head())
print("\nLMCs:", ranks['lmc'].unique())


import matplotlib.pyplot as plt
import seaborn as sns
from io import StringIO

# ── Already loaded ranks, now run GSEA per LMC ───────────────────────
lmc_list = sorted(ranks['lmc'].unique())

# ── Run prerank GSEA for each LMC ────────────────────────────────────
# Using Hallmarks gene sets — most interpretable
all_results = {}

for lmc in lmc_list:
    print(f'Running GSEA for {lmc}...')
    
    lmc_ranks = ranks[ranks['lmc'] == lmc][['gene', 'r']].copy()
    lmc_ranks = lmc_ranks.sort_values('r', ascending=False)
    lmc_ranks = lmc_ranks.drop_duplicates('gene')
    
    try:
        res = gp.prerank(
            rnk=lmc_ranks,
            gene_sets=[
                'MSigDB_Hallmark_2020',
                'Reactome_2022',
                'GO_Biological_Process_2023',
                'CellMarker_2024'
            ],
            outdir=None,
            min_size=10,
            max_size=500,
            permutation_num=1000,
            seed=42,
            verbose=False
        )
        all_results[lmc] = res.res2d
        n = len(res.res2d)
        print(f'  {lmc}: {n} gene sets tested')
    except Exception as e:
        print(f'  {lmc}: ERROR — {e}')

print('\nGSEA complete')

# ── Extract significant results ───────────────────────────────────────
sig_results = []
for lmc, df in all_results.items():
    df = df.copy()
    df['lmc'] = lmc
    sig = df[df['FDR q-val'] < 0.25]  # standard GSEA threshold
    sig_results.append(sig)
    print(f'{lmc}: {len(sig)} significant gene sets (FDR<0.25)')

sig_df = pd.concat(sig_results, ignore_index=True)
sig_df.to_csv('GSEA_Hallmarks_per_LMC.csv', index=False)
print(f'\nTotal significant results: {len(sig_df)}')
print(sig_df[['lmc', 'Term', 'NES', 'FDR q-val']].sort_values(['lmc', 'NES']).to_string())