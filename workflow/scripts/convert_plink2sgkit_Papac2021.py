from sgkit.io.plink import read_plink
import pandas as pd
import numpy as np
import xarray as xr

plink_prefix = snakemake.params['path']
file_meta_TS4 = snakemake.input['metadata_TS4']
file_meta_TS5 = snakemake.input['metadata_TS5']
file_meta_TS9 = snakemake.input['metadata_TS9']
file_name_mapping = snakemake.input['name_mapping']
file_anno = snakemake.input['anno']
bval_table = snakemake.input['bval_table']
rval_table = snakemake.input['rval_table']

zarr_store = snakemake.output['zarr_store']

variant_chunk_size = 10_000

cohorts_ref = [
	"Anatolia_Neolithic",
	"WHG",
	"Yamnaya",
]

ds = (
	read_plink(path=plink_prefix)
	.chunk({'alleles': 2, 'variants': variant_chunk_size})
)
ds = ds.isel(ploidy = [True, False])
ds = ds.isel(variants = (ds.variant_contig != 22) & (ds.variant_contig != 23)) # remove sex chr "23", "24"
ds.attrs['contigs'] = ds.attrs['contigs'][:-2]

TS4 = pd.read_csv(file_meta_TS4, sep='\t')
TS5 = pd.read_csv(file_meta_TS5, sep='\t')
name_mapping = pd.read_csv(file_name_mapping, sep='\t')
anno = pd.read_csv(file_anno, sep="\t")
anno = anno[anno['Version ID'].isin(ds.sample_id.values)]

samples = ds.sample_id.values


sample_group = [""]*ds.dims['samples']
for i, s in enumerate(samples):
	if np.isin(s, TS4['Sample']):
		w = np.where(s == TS4['Sample'])[0][0]
		sample_group[i] = TS4['qpAdm Grouping'].iloc[w]
	elif np.isin(s, TS5['Instance ID']):
		w = np.where(s == TS5['Instance ID'])[0][0]
		sample_group[i] = TS5['Group_ID'].iloc[w]

sample_group = [g if 'Yamnaya' not in g else 'Yamnaya' for g in sample_group ]
sample_group = [g if 'WHG' not in g else 'WHG' for g in sample_group ]
ds['sample_group'] = (['samples'], np.array(sample_group))


sample_date = np.full(ds.dims['samples'], np.nan, dtype=float)
for i, s in enumerate(samples):
	if np.isin(s, TS4['Sample']):
		w = np.where(s == TS4['Sample'])[0][0]
		sample_date[i] = TS4['BP (if directly dated)'].iloc[w]
		if np.isnan(sample_date[i]):
			estimate = TS4.at[w, '2-sigma low'] + (TS4.at[w, '2-sigma high'] - TS4.at[w, '2-sigma low'])/2
			sample_date[i] = estimate
	elif np.isin(s, anno['Version ID']):
		w = np.where(s == anno['Version ID'])[0][0]
		sample_date[i] = anno['Date mean in BP in years before 1950 CE [OxCal mu for a direct radiocarbon date, and average of range for a contextual date]'].iloc[w]
	else:
		print(f'sample nowhere to be found {s}')

ds['sample_date_bp'] = (['samples'], sample_date)


ds = ds.isel(samples=np.isin(samples, TS4['Sample']) | np.isin(ds.sample_group.values, cohorts_ref))


TS9 = pd.read_csv(file_meta_TS9, sep='\t')
TS9 = TS9.set_index('sample').join(name_mapping.set_index('library')).set_index('sample')
admix_values = np.array([
	list(TS9.loc[s, ["Anatolia_Neolithic", "WHG",  "Yamnaya_Samara"]])
	if np.isin(s, TS9.index) else [np.nan]*3
	for s in ds.sample_id.values
])
ds['sample_admixture'] = (
	['samples', 'cohorts_ref'],
	admix_values / 100
)
ds['cohorts_ref_id'] = (['cohorts_ref'], cohorts_ref)


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
