import sgkit
from sgkit.io.plink import read_plink
import pandas as pd
import numpy as np
import xarray as xr

plink_prefix = snakemake.params['path']
file_meta_S3 = snakemake.input['metadata_S3']
file_meta_S5 = snakemake.input['metadata_S5']
zarr_store = snakemake.output['zarr_store']

ds = read_plink(path=plink_prefix, chunks=5000)

meta_S3 = pd.read_csv(file_meta_S3, sep="\t")
meta_S5 = pd.read_csv(file_meta_S5, sep="\t")

meta = (
	pd.concat([meta_S3, meta_S5], join='outer', ignore_index=True)
	.drop_duplicates('Version ID', keep='last', ignore_index=True)
)

meta_sub = meta[np.isin(meta["Version ID"], ds.sample_id.values)].copy()

meta_sub.index = meta_sub['Version ID']
meta_sub = meta_sub.reindex(ds.sample_id.values)

filt1_col = 'Filter 1: Filter 0 plus p>0.01 qpAdm 3-way or qpAdm 2-way'
no_f1 = meta_sub[filt1_col].isna()
meta_sub.loc[no_f1, filt1_col] = meta_sub[no_f1]['Group ID']

ds['sample_group'] = xr.DataArray(meta_sub['Group ID'].to_numpy(dtype='str'), dims=['samples'])
ds['sample_date_bp'] = xr.DataArray(
	meta_sub['Date mean in BP [OxCal mu for a direct radiocarbon date, and average of range for a contextual date]'].to_numpy(),
	dims=['samples'],
)
ds['sample_snp_hits_autosomes'] = xr.DataArray(
	meta_sub['SNPs hit on autosomal targets'].to_numpy(),
	dims=['samples'],
)
ds['sample_filter_1'] = xr.DataArray(
	meta_sub[filt1_col].to_numpy(),
	dims=['samples'],
)
ds['sample_filter_0'] = xr.DataArray(
	meta_sub['Filter 0: >=30000 SNPs, PASS assessment, relevant period, not a first degree relative of higher coverage individuals)'].to_numpy(),
	dims=['samples'],
)
ds['sample_admixture'] = xr.DataArray(
	meta_sub.loc[:,['WHG (3-way model)', 'EEF (3-way model)', 'Steppe (3-way model)']].to_numpy(),
	dims=['samples', 'ancestries'],
)

sgkit.save_dataset(ds, store=zarr_store, mode='w')

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