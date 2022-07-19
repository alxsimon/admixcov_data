import sgkit
from sgkit.io.plink import read_plink
import pandas as pd
import numpy as np
import xarray as xr

plink_prefix = snakemake.params['path']
metadata = snakemake.input['metadata']
zarr_store = snakemake.output['zarr_store']

ds = read_plink(path=plink_prefix, chunks=5000)

meta = pd.read_csv(metadata, sep="\t")

meta_in = meta[np.isin(meta["Version ID"], ds.sample_id.values)].copy()
meta_out = meta[~np.isin(meta["Version ID"], ds.sample_id.values)].copy()

meta_in.index = meta_in['Version ID']
meta_in = meta_in.reindex(ds.sample_id.values)
ds['group_id'] = xr.DataArray(meta_in['Group ID'].to_numpy(dtype='str'), dims=['samples'])
ds['sample_date_bp'] = xr.DataArray(
	meta_in['Date mean in BP [OxCal mu for a direct radiocarbon date, and average of range for a contextual date]'].to_numpy(),
	dims=['samples'],
)
ds['sample_snp_hits_autosomes'] = xr.DataArray(
	meta_in['SNPs hit on autosomal targets'].to_numpy(),
	dims=['samples'],
)

sgkit.save_dataset(ds, store=zarr_store, mode='w')

meta_in.to_csv(
	snakemake.output['meta_in'],
	sep="\t",
	index=False,
)
meta_out.to_csv(
	snakemake.output['meta_out'],
	sep="\t",
	index=False,
)
