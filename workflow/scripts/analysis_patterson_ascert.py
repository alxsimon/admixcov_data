import sgkit as sg
import admixcov as ac
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import dask
import pickle
import pandas as pd


rng = np.random.default_rng(469829)

ds_full = sg.load_dataset(snakemake.input['zarr'])

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
ref_inds = np.isin(ds_full.sample_group.values, cohorts_ref)
kept_inds = (
        (ds_full.sample_filter_0.values == 'Use')
        & np.isin(ds_full.sample_filter_1.values, cohorts)
    ) | ref_inds
with dask.config.set(**{'array.slicing.split_large_chunks': True}):
    ds_full = ds_full.isel(samples=kept_inds)

ds_full['mask_cohorts'] = (
    ['cohorts', 'samples'],
    [(ds_full.sample_filter_1.values == p) for p in cohorts]
)
ds_full['mask_cohorts_ref'] = (
    ['cohorts_ref', 'samples'],
    [(ds_full.sample_group.values == p) for p in cohorts_ref]
)

ds_full['variant_count_nonmiss'] = (
    ['cohorts', 'variants'],
    np.stack([
        np.sum(~ds_full.call_genotype_mask.values[:, mask, 0], axis=1) 
        for mask in ds_full.mask_cohorts.values
    ])
)
ds_full['variant_count_nonmiss_ref'] = (
    ['cohorts_ref', 'variants'],
    np.stack([
        np.sum(~ds_full.call_genotype_mask.values[:, mask, 0], axis=1)
        for mask in ds_full.mask_cohorts_ref.values
    ])
)
kept_loci = (
    np.all((ds_full.variant_count_nonmiss.values > 10), axis=0) 
    & np.all((ds_full.variant_count_nonmiss_ref.values > 5), axis=0)
)
with dask.config.set(**{'array.slicing.split_large_chunks': False}):
    ds_full = ds_full.sel(variants=kept_loci)

ds_full['variant_chr_pos'] = (
    ['variants'],
    [str(chr + 1) + "_" + str(pos) for chr, pos in zip(
        ds_full.variant_contig.values,
        ds_full.variant_position.values
    )]
)


ascert = pd.concat((
        pd.read_csv(
            f"results/Patterson2022/ascertainment_info/Rohland2022_supp_chr{chr}.txt",
            sep='\t',
            header=None,
            usecols=[0, 1, 23],
            names=['chr', 'pos', 'ascertainment'],
        )
        for chr in range(1, 22)
    ),
    ignore_index=True,
)
ascert = ascert.loc[ascert.ascertainment.str.contains("HO=")]
ascert['HO'] = ascert.ascertainment.str.extract("=([01]+)")
ascert['chr_pos'] = ascert.chr.astype(str) + "_" + ascert.pos.astype(str)

A = np.array([[*map(int, list(str(x)))] for x in ascert.HO.to_numpy()], dtype=bool)

pop_ascert = {
    "1-French": 0,
    "2-Han": 1,
    "3-Papuan1": 2,
    "4-San": 3,
    "5-Yoruba": 4,
    "6-Mbuti": 5,
    "7-Karitiana": 6,
    "8-Sardinian": 7,
    "9-Melanesian": 8,
    "10-Cambodian": 9,
    "11-Mongolian": 10,
    "12-Papuan2": 11,
    "13-Denisova-San": 12,
}

focal_pops = [
    "1-French",
    "2-Han",
    "3-Papuan1",
    "4-San",
    "5-Yoruba",
    "13-Denisova-San",
]

results = []
for fp in focal_pops:
    tmp_ascert = ascert.loc[A[:, pop_ascert[fp]]]
    
    ds = ds_full.isel(
        variants=ds_full.variant_chr_pos.isin(tmp_ascert.chr_pos)
    )
    
    # need to update missing data
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
        G_nde.append(
            ac.get_matrix_sum(
                covmat[:i, :i] - admix_cov[:i, :i],
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
    straps_cov_nc = ac.bootstrap_stat(tiled_cov, weights, N_boot, rng=rng)
    straps_cov = ac.bootstrap_stat(tiled_corr_cov, weights, N_boot, rng=rng)

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
                rng=rng,
            )
        )
        straps_G.append(
            ac.bootstrap_ratio(
                np.stack([ac.get_matrix_sum(c) for c in tiled_corr_cov[:, :i, :i]]),
                tmp_totvar,
                weights,
                N_boot,
                rng=rng,
            )
        )
        straps_G_nc.append(
            ac.bootstrap_ratio(
                np.stack([ac.get_matrix_sum(c) for c in tiled_cov[:, :i, :i]]),
                tmp_totvar,
                weights,
                N_boot,
                rng=rng,
            )
        )
        straps_G_nde.append(
            ac.bootstrap_ratio(
                np.stack([ac.get_matrix_sum(c) for c in tiled_corr_cov_nde[:, :i, :i]]),
                tmp_totvar,
                weights,
                N_boot,
                rng=rng,
            )
        )
        straps_Ap.append(
            ac.bootstrap_ratio(
                np.stack([ac.get_matrix_sum(c, include_diag=True) for c in tiled_admix_cov[:, :i, :i]]),
                tmp_totvar,
                weights,
                N_boot,
                rng=rng,
            )
        )

    results.append({
        'times': np.array(times),
        'Q': Q,
        'straps_cov_nc': straps_cov_nc,
        'straps_cov': straps_cov,
        'straps_G': straps_G,
        'straps_G_nc': straps_G_nc,
        'straps_G_nde': straps_G_nde,
        'straps_Ap': straps_Ap,
        'N_snp': ds.dims['variants'],
    })


with open(snakemake.output['fig_data'], 'wb') as fw:
    pickle.dump(results, fw)
