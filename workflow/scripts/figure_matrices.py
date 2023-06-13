import admixcov as ac
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pickle

with open(snakemake.input['uk'], 'rb') as fr:
    uk_data = pickle.load(fr)

with open(snakemake.input['bo'], 'rb') as fr:
    bo_data = pickle.load(fr)

tri_lo = np.tril_indices(uk_data['straps_cov_nc'][1].shape[0], k=0)
for i in [0, 1, 2]: 
    uk_data['straps_cov_nc'][i][tri_lo] = np.nan

tri_lo = np.tril_indices(bo_data['straps_cov_nc'][1].shape[0], k=0)
for i in [0, 1, 2]: 
	bo_data['straps_cov_nc'][i][tri_lo] = np.nan

fig, axs = plt.subplots(1, 2, figsize=(10, 4))

labels = [f"$\\Delta p_{i}$" for i, t in enumerate(uk_data['times'][:-1])]
ac.plot_covmat_ci(uk_data['straps_cov_nc'], axs[0], delta_labels=labels)
axs[0].set_title("A", loc='left', fontdict={'fontweight': 'bold'})
axs[0].set_title('UK')

labels = [f"$\\Delta p_{i}$" for i, t in enumerate(bo_data['times'][:-1])]
ac.plot_covmat_ci(bo_data['straps_cov_nc'], axs[1], delta_labels=labels)
axs[1].set_title("B", loc='left', fontdict={'fontweight': 'bold'})
axs[1].set_title('Bohemia')

fig.tight_layout()
fig.savefig(snakemake.output['fig'])