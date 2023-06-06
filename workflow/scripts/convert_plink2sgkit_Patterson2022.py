from sgkit.io.plink import read_plink
import pandas as pd
import numpy as np
import xarray as xr

plink_prefix = snakemake.params['path']
file_meta_S3 = snakemake.input['metadata_S3']
file_meta_S5 = snakemake.input['metadata_S5']
file_modern_anno = snakemake.input['modern_anno']
bval_table = snakemake.input['bval_table']
rval_table = snakemake.input['rval_table']
zarr_store = snakemake.output['zarr_store']

variant_chunk_size = 10_000

cohorts = [
	'England.and.Wales_N',
	'England.and.Wales_C.EBA',
	'England.and.Wales_MBA',
	'England.and.Wales_LBA',
	'England.and.Wales_IA',
	'England.and.Wales_PostIA',
	'England.and.Wales_Modern',
]

cohorts_ref = [
	'WHGA',
	'Balkan_N',
	'OldSteppe',
]

ds = (
	read_plink(path=plink_prefix)
	.chunk({'alleles': 2, 'variants': variant_chunk_size})
)

ds['cohorts_id'] = (['cohorts'], cohorts)
ds['cohorts_ref_id'] = (['cohorts_ref'], cohorts_ref)

ds = ds.isel(ploidy = [True, False])

ds = ds.isel(variants = (ds.variant_contig != 22) & (ds.variant_contig != 23)) # remove sex chr "23", "24"
ds.attrs['contigs'] = ds.attrs['contigs'][:-2]

# integrate metadata
meta_S3 = pd.read_csv(file_meta_S3, sep="\t")
meta_S5 = pd.read_csv(file_meta_S5, sep="\t")
modern_anno = pd.read_csv(file_modern_anno, sep="\t")

modern_anno.rename(
	columns={'Date mean in BP in years before 1950 CE [OxCal mu for a direct radiocarbon date, and average of range for a contextual date]': 'Date mean in BP [OxCal mu for a direct radiocarbon date, and average of range for a contextual date]'},
	inplace=True,
) # correct one important column name to fit other metadata

meta = (
	pd.concat([modern_anno, meta_S3, meta_S5], join='outer', ignore_index=True)
	.drop_duplicates('Version ID', keep='last', ignore_index=True)
)

meta_sub = meta[np.isin(meta["Version ID"], ds.sample_id.values)].copy()

meta_sub.index = meta_sub['Version ID']
meta_sub = meta_sub.reindex(ds.sample_id.values)

filt1_col = 'Filter 1: Filter 0 plus p>0.01 qpAdm 3-way or qpAdm 2-way'
no_f1 = meta_sub[filt1_col].isna()
meta_sub.loc[no_f1, filt1_col] = meta_sub[no_f1]['Group ID']

ds['sample_group'] = (['samples'], meta_sub['Group ID'].to_numpy(dtype='str'))

ds['sample_date_bp'] = (
	['samples'],
	meta_sub['Date mean in BP [OxCal mu for a direct radiocarbon date, and average of range for a contextual date]'].to_numpy(),
)
ds['sample_snp_hits_autosomes'] = (
	['samples'],
	meta_sub['SNPs hit on autosomal targets'].to_numpy(),
)
ds['sample_filter_1'] = (
	['samples'],
	meta_sub[filt1_col].to_numpy(),
)
ds['sample_filter_0'] = (
	['samples'],
	meta_sub['Filter 0: >=30000 SNPs, PASS assessment, relevant period, not a first degree relative of higher coverage individuals)'].to_numpy(),
)
ds['sample_admixture'] = (
	['samples', 'cohorts_ref'],
	meta_sub.loc[:,['WHG (3-way model)', 'EEF (3-way model)', 'Steppe (3-way model)']].to_numpy(),
)

sample_cohort = np.array([
	np.where(np.array(cohorts) == x)[0][0] if np.isin(x, cohorts) & (y == 'Use') else -9
	for x, y in zip(ds.sample_filter_1.values, ds.sample_filter_0.values)
], dtype=int)
ds['sample_cohort'] = (['samples'], sample_cohort)

# integrate b-value and recombination rate
ds = ds.unify_chunks()

bmap = pd.read_csv(bval_table, sep="\t")
bmap['bval'] = bmap['bval'] * 1e-3
rmap = pd.read_csv(rval_table, sep="\t")

def get_val(df, v_name, chr_idx, pos):
	val = df[v_name][(df.chr_num.values == chr_idx) & (df.start.values <= pos) & (df.end.values > pos)].values
	if val.size == 0:
		return np.nan
	else:
		return val[0]

def get_value_from_intervals(ds, df, v_name):
	df['chr_num'] = (df.chr.str.replace("chr", "").values.astype(int) - 1)
	snp_val = [
		get_val(df, v_name, chr_idx, pos)
		for chr_idx, pos in zip(ds.variant_contig.values, ds.variant_position.values)
	]
	return xr.DataArray(data=snp_val, dims=('variants'), name=f'variant_{v_name}')

da_bval = xr.map_blocks(get_value_from_intervals, ds, kwargs={'df': bmap, 'v_name': 'bval'}).compute()
da_rmap = xr.map_blocks(get_value_from_intervals, ds, kwargs={'df': rmap, 'v_name': 'rate'}).compute()

ds['variant_bval'] = da_bval
ds['variant_rate'] = da_rmap

# save
ds = ds.chunk(variant_chunk_size)

for var in ds:
	ds[var].encoding.clear()
for var in ds:
	if ds[var].dtype == 'object':
		max_len_names = max([len(x) if x is not np.nan else 0 for x in ds[var].values])
		ds[var] = ds[var].astype(f"U{max_len_names}")

ds.to_zarr(store=zarr_store, mode='w')

meta_sub.to_csv(
	snakemake.output['meta_out'],
	sep="\t",
	index=False,
)

# Create ind file useful for smartpca
ind = pd.DataFrame(meta_sub['Version ID'])
ind['sex'] = 'U'
ind['group'] = meta_sub[filt1_col]
ind.to_csv(snakemake.output['out_ind'], sep="\t", index=False, header=False)