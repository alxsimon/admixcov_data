import sgkit as sg
import admixcov as ac
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import dask

ds = sg.load_dataset(snakemake.input['zarr'])

# grouping
time_periods_grouping = {
	'Bohemia_N': ['Bohemia_N'],
	'Bohemia_PE': ['Bohemia_PE'],
	'Bohemia_EE': ['Bohemia_EE'],
	'Bohemia_ME': ['Bohemia_ME'], # ['Bohemia_ME_Baden', 'Bohemia_ME_Rivnac', 'Bohemia_ME_GAC'],
	'Bohemia_CW': ['Bohemia_CW_Early', 'Bohemia_CW_Late'],
	'Bohemia_BB': ['Bohemia_BB_Early', 'Bohemia_BB_Late'],
	'Bohemia_Unetice': ['Bohemia_Unetice_preClassical', 'Bohemia_Unetice_Classical'],
}

cohorts = list(time_periods_grouping.keys())
cohorts_ref = [
	"Anatolia_Neolithic",
	"WHG",
	"Yamnaya",
]

ds['mask_cohorts'] = (
	['cohorts', 'samples'],
	[[i in gs for i in ds.sample_group.values] for gs in time_periods_grouping.values()]
)
ds['mask_cohorts_ref'] = (
	['cohorts_ref', 'samples'],
	[(ds.sample_group.values == p) for p in cohorts_ref]
)

ds['variant_count_nonmiss'] = (
	['cohorts', 'variants'],
	np.stack([
		np.sum(~ds.call_genotype_mask.values[:, mask, 0], axis=1) 
		for mask in ds.mask_cohorts.values
	])
)
ds['variant_count_nonmiss_ref'] = (
	['cohorts_ref', 'variants'],
	np.stack([
		np.sum(~ds.call_genotype_mask.values[:, mask, 0], axis=1)
		for mask in ds.mask_cohorts_ref.values
	])
)
kept_loci = (
	np.all((ds.variant_count_nonmiss.values > 2), axis=0) 
	& np.all((ds.variant_count_nonmiss_ref.values > 2), axis=0)
)
with dask.config.set(**{'array.slicing.split_large_chunks': False}):
	ds = ds.sel(variants=kept_loci)

times = [np.mean(ds.sample_date_bp.values[mask]) for mask in ds.mask_cohorts.values]

geno = ds.call_genotype.values[:,:,0].T.astype(float)
geno[geno == -1] = np.nan

ref_af = np.stack([np.nanmean(geno[mask], axis=0) for mask in ds.mask_cohorts_ref.values])
sample_size_ref = [np.sum(mask) for mask in ds.mask_cohorts_ref.values]

def allele_freq(geno, mask):
	return np.nanmean(geno[mask], axis=0)

af = np.stack([
	allele_freq(geno, mask)
	for mask in ds.mask_cohorts.values
])

sample_size = np.array([
	np.sum(mask)
	for mask in ds.mask_cohorts.values
])

Q = admix = np.stack([
	np.mean(ds.sample_admixture[mask], axis=0)
	for mask in ds.mask_cohorts.values
])

covmat = ac.get_covariance_matrix(
	af,
	bias=True,
	sample_size=ds.variant_count_nonmiss.values,
)
admix_cov = ac.get_admix_covariance_matrix(
	Q,
	ref_af=ref_af,
	bias=True,
	ref_sample_size=ds.variant_count_nonmiss_ref.values,
)

alpha_mask = np.array([ # Anatolia, WHG, Yamnaya
	[0, 1, 0],
	[0, 0, 0],
	[0, 0, 0],
	[0, 0, 1],
	[1, 0, 0],
	[0, 0, 1],
], dtype=bool) # says which alpha is different from zero
alphas = ac.q2a_simple(Q, alpha_mask)
var_drift = ac.solve_for_variances(
	np.diag(covmat - admix_cov),
	alphas,
)
drift_err = ac.get_drift_err_matrix(var_drift, alphas)

k = covmat.shape[0]
G = []
G_nc = []
Ap = []
totvar = []
for i in range(1, k + 1):
    totvar.append(np.sum(covmat[:i, :i]))
    G.append(
        ac.get_matrix_sum(
            covmat[:i, :i] - admix_cov[:i, :i] - drift_err[:i, :i],
            include_diag=False, abs=False
        ) / totvar[-1]
    )
    G_nc.append(
        ac.get_matrix_sum(
            covmat[:i, :i],
            include_diag=False, abs=False
        ) / totvar[-1]
    )
    Ap.append(
        ac.get_matrix_sum(
            admix_cov[:i, :i],
            include_diag=True, abs=False
        ) / totvar[-1]
    )
    
N_boot = 1e4
tile_idxs = ac.sg.create_tile_idxs(ds, type='variant', size=1500)
sizes = [x.size for x in tile_idxs] # Number of SNPs in tiles 

n_sample = ds.variant_count_nonmiss.values
n_sample_ref = ds.variant_count_nonmiss_ref.values
tiled_af = [af[:, idx] for idx in tile_idxs]
tiled_sample_size = [n_sample[:, idx] for idx in tile_idxs]

assert af.shape == n_sample.shape
tiled_cov = np.stack([
    ac.get_covariance_matrix(a, bias=True, sample_size=n)
    for a, n in zip(tiled_af, tiled_sample_size)
])

assert ref_af.shape == n_sample_ref.shape
tiled_admix_cov = np.stack([
    ac.get_admix_covariance_matrix(
        Q, ref_af[:, idx], bias=True,
        ref_sample_size=n_sample_ref[:, idx],
    )
    for idx in tile_idxs
])

tiled_drift_err = [
    ac.get_drift_err_matrix(
        ac.solve_for_variances(np.diag(c - a), alphas),
        alphas,
    )
    for c, a in zip(tiled_cov, tiled_admix_cov)
]

tiled_corr_cov = np.stack([
    c - a - d for c, a, d in zip(tiled_cov, tiled_admix_cov, tiled_drift_err)
])

n_loci = np.array([tile.size for tile in tile_idxs])
weights = n_loci / np.sum(n_loci)

# do the bootstraps
straps_cov_nc = ac.bootstrap_stat(tiled_cov, weights, N_boot)
straps_cov = ac.bootstrap_stat(tiled_corr_cov, weights, N_boot)

straps_G = []
straps_G_nc = []
straps_Ap = []
straps_totvar = []
k = tiled_cov.shape[1]
for i in range(1, k + 1):
    tmp_totvar = np.sum(tiled_cov[:, :i, :i], axis=(1, 2))
    straps_totvar.append(
        ac.bootstrap_stat(
            tmp_totvar,
            weights,
            N_boot,
        )
    )
    straps_G.append(
        ac.bootstrap_ratio(
            np.stack([ac.get_matrix_sum(c) for c in tiled_corr_cov[:, :i, :i]]),
            tmp_totvar,
            weights,
            N_boot,
            statistic=G[i - 1],
        )
    )
    straps_G_nc.append(
        ac.bootstrap_ratio(
            np.stack([ac.get_matrix_sum(c) for c in tiled_cov[:, :i, :i]]),
            tmp_totvar,
            weights,
            N_boot,
            statistic=G_nc[i - 1],
        )
    )
    straps_Ap.append(
        ac.bootstrap_ratio(
            np.stack([ac.get_matrix_sum(c, include_diag=True) for c in tiled_admix_cov[:, :i, :i]]),
            tmp_totvar,
            weights,
            N_boot,
            statistic=Ap[i - 1],
        )
    )


# Plotting

time_padding = 100

colors_oi = [
    '#000000', # black
    '#D55E00', # vermillion
    '#0072B2', # blue
    '#009E73', # green
    '#E69F00', # orange
    '#56B4E9', # sky blue
    '#CC79A7', # pink
    '#F0E442', # yellow
]

times = np.array(times) # ensure it is an array

fig, axs = plt.subplots(3, 2, figsize=(10, 8))

fmts = ['-o', '-s', '-^']
for i, pop in enumerate(ds.cohorts_ref_id.values):
	axs[0,0].plot(times, Q[:,i], fmts[i], label=pop, color=colors_oi[i])
axs[0,0].set_xlim(times[0] + time_padding, times[-1] - time_padding)
axs[0,0].set_ylim(top=1)
axs[0,0].set_ylabel("Mean ancestry")
axs[0,0].set_xlabel("Time point")
axs[0,0].legend(loc='lower center', bbox_to_anchor=(0.5, 1), ncol=3)

combined_ci = ac.combine_covmat_CIs(straps_cov, straps_cov_nc)
scale_max = np.max(np.abs([np.nanmin(combined_ci[1] - np.diag(np.diag(combined_ci[1]))), np.nanmax(combined_ci[1] - np.diag(np.diag(combined_ci[1])))]))
ac.plot_covmat_ci(combined_ci, axs[0, 1], scale_max)
axs[0,1].set_title('covariance matrix (raw lower, corrected upper)')

ac.cov_lineplot(times, straps_cov_nc, axs[1, 0], colors=colors_oi, time_padding=time_padding, d=50, marker='o')
axs[1, 0].set_ylabel('raw covariance (without bias)')
ac.cov_lineplot(times, straps_cov, axs[1, 1], colors=colors_oi, time_padding=time_padding, d=50, marker='o', ylim=axs[1, 0].get_ylim())
axs[1, 1].set_ylabel('admixture corrected covariance')

ac.plot_ci_line(times[1:], np.stack(straps_totvar).T, ax=axs[2, 0], color='black', marker='o')
axs[2, 0].set_xlim(times[1] + time_padding, times[-1] - time_padding)
axs[2, 0].set_ylim(0)
axs[2, 0].set_ylabel('Total variance (t)')

x_shift = 50
ac.plot_ci_line(times[1:] + x_shift, np.stack(straps_G_nc).T, ax=axs[2, 1], linestyle='dashed', marker='o')
ac.plot_ci_line(times[1:], np.stack(straps_G).T, ax=axs[2, 1], marker='o')
ac.plot_ci_line(times[1:] - x_shift, np.stack(straps_Ap).T, ax=axs[2, 1], color='blue', marker='s')
axs[2, 1].set_xlim(times[1] + time_padding, times[-1] - time_padding)
axs[2, 1].hlines(y=0, xmin=times[-1] - time_padding, xmax=times[1] + time_padding, colors='black', linestyles='dotted')
axs[2, 1].set_xlabel('time')
axs[2, 1].set_ylabel("G(t) or A'(t)")

fig.tight_layout()

fig.savefig(snakemake.output['fig'])