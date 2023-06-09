import sgkit as sg
import admixcov as ac
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import dask
import pickle

report = open(snakemake.output['report'], 'w')

ds = sg.load_dataset(snakemake.input['zarr'])
report.write("Dataset import:\n")
print(ds.dims, file=report)

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
ref_inds = np.isin(ds.sample_group.values, cohorts_ref)
kept_inds = (
        (ds.sample_filter_0.values == 'Use')
        & np.isin(ds.sample_filter_1.values, cohorts)
    ) | ref_inds
with dask.config.set(**{'array.slicing.split_large_chunks': True}):
    ds = ds.isel(samples=kept_inds)

ds['mask_cohorts'] = (
    ['cohorts', 'samples'],
    [(ds.sample_filter_1.values == p) for p in cohorts]
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
    np.all((ds.variant_count_nonmiss.values > 10), axis=0) 
    & np.all((ds.variant_count_nonmiss_ref.values > 5), axis=0)
)
with dask.config.set(**{'array.slicing.split_large_chunks': False}):
    ds = ds.sel(variants=kept_loci)
report.write("\n==========\nFiltered dataset:\n")
print(ds.dims, file=report)

times = [np.mean(ds.sample_date_bp.values[mask]) for mask in ds.mask_cohorts.values]
report.write("\n==========\nMean times:\n")
print(times, file=report)

report.write("\n==========\nSample sizes:\n")
print(f"Samples: {[x.sum() for x in ds.mask_cohorts.values]}", file=report)
print(f"for {cohorts}", file=report)
print(f"Refs: {[x.sum() for x in ds.mask_cohorts_ref.values]}", file=report)
print(f"for {cohorts_ref}", file=report)

report.write("\n==========\nNon-missing data means:\n")
print(f"Samples: {np.mean(ds.variant_count_nonmiss.values, axis=1)}", file=report)
print(f"Refs: {np.mean(ds.variant_count_nonmiss_ref.values, axis=1)}", file=report)

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

Q  = np.stack([
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

alpha_mask = np.array([ # WHG, EEF, Steppe
    [0, 0, 1],
    [0, 1, 0],
    [0, 1, 0],
    [0, 1, 0],
    [1, 0, 0],
    [0, 1, 0],
], dtype=bool) # says which alpha is different from zero
alphas = ac.q2a_simple(Q, alpha_mask)
var_drift = ac.solve_for_variances(
    np.diag(covmat - admix_cov),
    alphas,
)
drift_err = ac.get_drift_err_matrix(var_drift, alphas)

report.write("\n==========\nCovariances and drift component:\n")
report.write("covmat:\n")
print(covmat, file=report)
report.write("admix_cov:\n")
print(admix_cov, file=report)
report.write("Drift err:\n")
print(drift_err, file=report)


k = covmat.shape[0]
G = []
G_nc = []
G_nde = []
Ap = []
totvar = []
for i in range(1, k + 1):
    totvar.append(np.sum(covmat[:i, :i]))
    G.append(
        ac.get_matrix_sum(
            covmat[:i, :i] - admix_cov[:i, :i],
            include_diag=False, abs=False
        ) / totvar[-1]
    )
    G_nc.append(
        ac.get_matrix_sum(
            covmat[:i, :i],
            include_diag=False, abs=False
        ) / totvar[-1]
    )
    G_nde.append(
        ac.get_matrix_sum(
            covmat[:i, :i] - admix_cov[:i, :i] - drift_err[:i, :i],
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
tile_idxs = ac.sg.create_tile_idxs(ds, type='variant', size=1_000)
# tile_idxs = ac.sg.create_tile_idxs(ds, type='position', size=1_000_000)
# tile_idxs = [t for t in tile_idxs if len(t) >= 100] # filter
sizes = [x.size for x in tile_idxs] # Number of SNPs in tiles

report.write("\n==========\nBootstrap info:\n")
report.write("number of tiles:\n")
print(len(tile_idxs), file=report)
report.write("min, max, median sizes:\n")
print(f"{np.min(sizes)}, {np.max(sizes)}, {np.median(sizes)}", file=report)

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

tiled_drift_err = np.stack([
    ac.get_drift_err_matrix(
        ac.solve_for_variances(np.diag(c - a), alphas),
        alphas,
    )
    for c, a in zip(tiled_cov, tiled_admix_cov)
])

tiled_corr_cov_nde = np.stack([
    c - a for c, a in zip(tiled_cov, tiled_admix_cov)
])
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
straps_G_nde = []
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
        )
    )
    straps_G_nc.append(
        ac.bootstrap_ratio(
            np.stack([ac.get_matrix_sum(c) for c in tiled_cov[:, :i, :i]]),
            tmp_totvar,
            weights,
            N_boot,
        )
    )
    straps_G_nde.append(
        ac.bootstrap_ratio(
            np.stack([ac.get_matrix_sum(c) for c in tiled_corr_cov_nde[:, :i, :i]]),
            tmp_totvar,
            weights,
            N_boot,
        )
    )
    straps_Ap.append(
        ac.bootstrap_ratio(
            np.stack([ac.get_matrix_sum(c, include_diag=True) for c in tiled_admix_cov[:, :i, :i]]),
            tmp_totvar,
            weights,
            N_boot,
        )
    )

report.write("\n==========\nConfidence intervals:\n")
report.write("G:\n")
print(straps_G, file=report)
report.write("G_nc:\n")
print(straps_G_nc, file=report)
report.write("G_nde:\n")
print(straps_G_nde, file=report)
report.write("Ap:\n")
print(straps_Ap, file=report)


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

fig, axs = plt.subplots(2, 2, figsize=(10, 8), layout='constrained')

k, l = (0, 1)
fmts = ['-o', '-s', '-^']
labels = ['WHG-like', 'EEF-like', 'Steppe-like']
for i, pop in enumerate(ds.cohorts_ref_id.values):
    axs[k, l].plot(times, Q[:,i], fmts[i], label=labels[i], color=colors_oi[i])
axs[k, l].set_xlim(times[0] + time_padding, times[-1] - time_padding)
axs[k, l].set_ylim(top=1)
axs[k, l].set_ylabel("Mean ancestry")
axs[k, l].set_xlabel("Time (years BP)")
axs[k, l].legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=3)
axs[k, l].set_title("B", loc='left', fontdict={'fontweight': 'bold'})

# combined_ci = ac.combine_covmat_CIs(straps_cov, straps_cov_nc)
# scale_max = np.max(np.abs([np.nanmin(combined_ci[1] - np.diag(np.diag(combined_ci[1]))), np.nanmax(combined_ci[1] - np.diag(np.diag(combined_ci[1])))]))
# ac.plot_covmat_ci(combined_ci, axs[0, 1], scale_max)
# axs[0,1].set_title('covariance matrix (raw lower, corrected upper)')

x_shift = 50
k, l = (0, 0)
ac.cov_lineplot(times, straps_cov_nc, axs[k, l], colors=colors_oi, time_padding=time_padding, d=x_shift, marker='o')
axs[k, l].set_xlim(times[1] + x_shift + time_padding, times[-2] - x_shift - time_padding)
axs[k, l].set_ylabel("Cov($\\Delta p_i$, $\\Delta p_t$)")
axs[k, l].set_xlabel("t")
axs[k, l].set_title('Before admix. correction')
axs[k, l].set_title("A", loc='left', fontdict={'fontweight': 'bold'})

k, l = (1, 0)
ac.cov_lineplot(times, straps_cov, axs[k, l], colors=colors_oi, time_padding=time_padding, d=x_shift, marker='o', ylim=axs[0, 0].get_ylim())
axs[k, l].set_xlim(times[1] + x_shift + time_padding, times[-2] - x_shift - time_padding)
axs[k, l].set_ylabel("Cov($\\Delta p_i$, $\\Delta p_t$)")
axs[k, l].set_xlabel('t')
axs[k, l].set_title('After admix. correction')
axs[k, l].set_title("C", loc='left', fontdict={'fontweight': 'bold'})
axs[k, l].legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), title="$\\Delta p_i$", ncol=3)

# ac.plot_ci_line(times[1:], np.stack(straps_totvar).T, ax=axs[k, l], color='black', marker='o')
# axs[k, l].set_xlim(times[1] + time_padding, times[-1] - time_padding)
# axs[k, l].set_ylim(0)
# axs[k, l].set_ylabel('Total variance (t)')

k, l = (1, 1)
ac.plot_ci_line(times[1:] + x_shift, np.stack(straps_G_nc).T, ax=axs[k, l], linestyle='dashed', marker='o', label='$G_{nc}$')
ac.plot_ci_line(times[1:] + 2 * x_shift, np.stack(straps_G_nde).T, ax=axs[k, l], linestyle='dashdot', marker='^', label='$G_{nde}$')
ac.plot_ci_line(times[1:], np.stack(straps_G).T, ax=axs[k, l], marker='o', label='G')
ac.plot_ci_line(times[1:] - x_shift, np.stack(straps_Ap).T, ax=axs[k, l], color='blue', marker='s', label='A\'')
axs[k, l].set_xlim(times[1] + x_shift + time_padding, times[-1] - x_shift - time_padding)
axs[k, l].hlines(y=0, xmin=times[-1] - time_padding, xmax=times[1] + time_padding, colors='grey', linestyles='dotted')
axs[k, l].set_xlabel('t')
axs[k, l].set_ylabel("Proportion of variance ($p_t - p_{t0}$)")
axs[k, l].legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=4)
axs[k, l].set_title("D", loc='left', fontdict={'fontweight': 'bold'})
for ci, t in zip(straps_G, times[1:]):
    if ci[0]*ci[2] > 0:
        axs[k, l].annotate("*", xy=(t, 0.1))

fig.savefig(snakemake.output['fig'])


#==========
# Do some reduced bootstrap to have a genome wide distribution
# to compare to

rng = np.random.default_rng()
tiled_corr_totvar = np.sum(tiled_corr_cov, axis=(1, 2))

# do the bootstraps
L = len(tiled_corr_totvar)
subset_L = round(len(tile_idxs)/5)
straps = []
for _ in np.arange(N_boot):
    bidx = rng.integers(0, L, size=subset_L)
    straps.append(
        ac.weighted_mean(
            tiled_corr_totvar[bidx],
            weights[bidx],
            axis=0,
        )
    )
straps = np.array(straps)
That = np.mean(straps, axis=0)
totvar_corr_CI_sub = ac.bootstrap_ci(That, straps, alpha=0.05, axis=0)


#================
# Binning
# recombination
nanmask = ~np.isnan(ds.variant_rate.values)
ds_rate = ds.sel(variants=nanmask)

q = np.quantile(ds_rate.variant_rate.values, list(np.arange(1/5, 1, 1/5)))
report.write("\n==========\nRecombination quantiles:\n")
min, max = (
    np.min(ds_rate.variant_rate.values),
    np.max(ds_rate.variant_rate.values)
)
print(f"min: {min}", file=report)
print(f"max: {max}", file=report)
print(q, file=report)

bins = np.digitize(ds_rate.variant_rate.values, q)

bin_res = []
for bin in np.unique(bins):
	bin_res.append(
		ac.sg.ds2stats(
			ds_rate.sel(variants=(bins == bin)),
			alpha_mask,
            tile_type='variant',
			tile_size=1_000,
		)
	)

# 
G_CI = [x[5] for x in bin_res]
Ap_CI = [x[6] for x in bin_res]

fig1, axs1 = plt.subplots(1, 2, figsize=(10, 5), layout='constrained') # G
fig2, axs2 = plt.subplots(1, 2, figsize=(10, 5), layout='constrained') # individual Var
fig3, axs3 = plt.subplots(2, 2, figsize=(10, 10), layout='constrained') # totvar

ac.plot_ci_line(np.unique(bins), np.stack(G_CI).T, axs1[0], marker='o', label='G')
ac.plot_ci_line(np.unique(bins), np.stack(Ap_CI).T, axs1[0], marker='o', color='b', label='A\'')
axs1[0].hlines(y=0, xmin=0, xmax=4, color='black', linestyles='dotted')
axs1[0].set_xlabel('Recombination bin')
axs1[0].set_ylabel('Proportion of variance')

covmats = [x[4] for x in bin_res]
# variances devided by half hz
vardiag_bins_CIs = [
	(
		np.array([C[0][i,i] / bin_res[j][13][i] for j, C in enumerate(covmats)]),
		np.array([C[1][i,i] / bin_res[j][13][i] for j, C in enumerate(covmats)]),
		np.array([C[2][i,i] / bin_res[j][13][i] for j, C in enumerate(covmats)]),
	)
	for i in range(6)]
for i, ci in enumerate(vardiag_bins_CIs):
	ac.plot_ci_line(np.unique(bins) + 0.1 * i, ci, axs2[0], marker='o', label=f'p_{int(times[i])}', color=colors_oi[i])
axs2[0].set_xlabel('Recombination bin')
axs2[0].set_ylabel('$Var(\Delta p_t) / p_t(1 - p_t)$')
axs2[0].hlines(y=0, xmin=0, xmax=4, color='black', linestyles='dotted')

#   uncorrected
# totvar
ac.plot_ci_line(
    np.unique(bins),
    np.stack([x[7] for x in bin_res]).T,
    axs3[0, 0], marker='o',
    color=colors_oi[0],
    label='Total variance',
)
axs3[0, 0].set_xlabel('Recombination bin')
axs3[0, 0].xaxis.set_major_locator(ticker.MultipleLocator(1))

# sum var
ac.plot_ci_line(
    np.unique(bins) + 0.1,
    np.stack([x[9] for x in bin_res]).T,
    axs3[0, 0], marker='s',
    color=colors_oi[1],
    label='Sum of variances',
)

# sum cov
ac.plot_ci_line(
    np.unique(bins) + 0.2,
    np.stack([x[11] for x in bin_res]).T,
    axs3[0, 0], marker='^',
    color=colors_oi[2],
    label='Sum of covariances',
)

#   corrected
# totvar
ac.plot_ci_line(
    np.unique(bins),
    np.stack([x[8] for x in bin_res]).T,
    axs3[1, 0], marker='o',
    color=colors_oi[0],
    label='Total variance',
)
axs3[1, 0].set_xlabel('Recombination bin')
axs3[1, 0].xaxis.set_major_locator(ticker.MultipleLocator(1))

# sum var
ac.plot_ci_line(
    np.unique(bins) + 0.1,
    np.stack([x[10] for x in bin_res]).T,
    axs3[1, 0], marker='s',
    color=colors_oi[1],
    label='Sum of variances',
)

# sum cov
ac.plot_ci_line(
    np.unique(bins) + 0.2,
    np.stack([x[12] for x in bin_res]).T,
    axs3[1, 0], marker='^',
    color=colors_oi[2],
    label='Sum of covariances',
)

#==========
# bval
nanmask = ~np.isnan(ds.variant_bval.values)
ds_bval = ds.sel(variants=nanmask)

q = np.quantile(ds_bval.variant_bval.values, list(np.arange(1/5, 1, 1/5)))
report.write("\n==========\nB-values quantiles:\n")
min, max = (
    np.min(ds_rate.variant_bval.values),
    np.max(ds_rate.variant_bval.values)
)
print(f"min: {min}", file=report)
print(f"max: {max}", file=report)
print(q, file=report)

bins = np.digitize(ds_bval.variant_bval.values, q)

bin_res = []
for bin in np.unique(bins):
	bin_res.append(
		ac.sg.ds2stats(
			ds_bval.sel(variants=(bins == bin)),
			alpha_mask,
            tile_type='variant',
			tile_size=1_000,
		)
	)

# 
G_CI = [x[5] for x in bin_res]
Ap_CI = [x[6] for x in bin_res]

ac.plot_ci_line(np.unique(bins), np.stack(G_CI).T, axs1[1], marker='o', label='G')
ac.plot_ci_line(np.unique(bins), np.stack(Ap_CI).T, axs1[1], marker='o', color='b', label='A\'')
axs1[1].hlines(y=0, xmin=0, xmax=4, color='black', linestyles='dotted')
axs1[1].set_xlabel('B-value bin')
axs1[1].set_ylabel('Proportion of variance')

covmats = [x[4] for x in bin_res]
# variances devided by half hz
vardiag_bins_CIs = [
	(
		np.array([C[0][i,i] / bin_res[j][13][i] for j, C in enumerate(covmats)]),
		np.array([C[1][i,i] / bin_res[j][13][i] for j, C in enumerate(covmats)]),
		np.array([C[2][i,i] / bin_res[j][13][i] for j, C in enumerate(covmats)]),
	)
	for i in range(len(times) - 1)]
for i, ci in enumerate(vardiag_bins_CIs):
	ac.plot_ci_line(np.unique(bins) + 0.1 * i, ci, axs2[1], marker='o', label=f'p_{int(times[i])}', color=colors_oi[i])
axs2[1].set_xlabel('B-value bin')
axs2[1].set_ylabel('$Var(\Delta p_t) / p_t(1 - p_t)$')
axs2[1].hlines(y=0, xmin=0, xmax=4, color='black', linestyles='dotted')

#   uncorrected
# totvar
ac.plot_ci_line(
    np.unique(bins),
    np.stack([x[7] for x in bin_res]).T,
    axs3[0, 1], marker='o',
    color=colors_oi[0],
    label='Total variance',
)
axs3[0, 1].set_xlabel('B-value bin')
axs3[0, 1].xaxis.set_major_locator(ticker.MultipleLocator(1))

# sum var
ac.plot_ci_line(
    np.unique(bins) + 0.1,
    np.stack([x[9] for x in bin_res]).T,
    axs3[0, 1], marker='s',
    color=colors_oi[1],
    label='Sum of variances',
)

# sum cov
ac.plot_ci_line(
    np.unique(bins) + 0.2,
    np.stack([x[11] for x in bin_res]).T,
    axs3[0, 1], marker='^',
    color=colors_oi[2],
    label='Sum of covariances',
)

#   corrected
# totvar

axs3[1, 1].hlines(y=totvar_corr_CI_sub[1], xmin=0, xmax=4, linestyles='dashdot')
axs3[1, 1].fill_between(x=[0, 4], y1=[totvar_corr_CI_sub[0]]*2, y2=[totvar_corr_CI_sub[2]]*2, color='b', alpha=.15)

ac.plot_ci_line(
    np.unique(bins),
    np.stack([x[8] for x in bin_res]).T,
    axs3[1, 1], marker='o',
    color=colors_oi[0],
    label='Total variance',
)
axs3[1, 1].set_xlabel('B-value bin')
axs3[1, 1].xaxis.set_major_locator(ticker.MultipleLocator(1))

# sum var
ac.plot_ci_line(
    np.unique(bins) + 0.1,
    np.stack([x[10] for x in bin_res]).T,
    axs3[1, 1], marker='s',
    color=colors_oi[1],
    label='Sum of variances',
)

# sum cov
ac.plot_ci_line(
    np.unique(bins) + 0.2,
    np.stack([x[12] for x in bin_res]).T,
    axs3[1, 1], marker='^',
    color=colors_oi[2],
    label='Sum of covariances',
)

#==============

handles, labels = axs1[0].get_legend_handles_labels()
fig1.legend(handles, labels, loc='outside upper center', ncols=2)
fig1.savefig(snakemake.output['fig_bins_G'])

handles, labels = axs2[0].get_legend_handles_labels()
fig2.legend(handles, labels, loc='outside right upper', title='$p_t$')
fig2.savefig(snakemake.output['fig_bins_var'])

handles, labels = axs3[0, 0].get_legend_handles_labels()
fig3.legend(handles, labels, loc='outside upper center')
fig3.savefig(snakemake.output['fig_bins_totvar'])

#==============
matrix_data = {
    'straps_cov_nc': straps_cov_nc,
    'times': times,
}
with open(snakemake.output['matrix_data'], 'wb') as fw:
    pickle.dump(matrix_data, fw)


#==============
# close report
report.close()