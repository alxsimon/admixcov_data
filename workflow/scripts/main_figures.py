import admixcov as ac
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import pickle

import matplotlib.ticker as tkr
loc = tkr.MultipleLocator(base=1.0)
# sci notation formatter
formatter = tkr.ScalarFormatter(useMathText=True)
formatter.set_scientific(True)
formatter.set_powerlimits((0, 0))

with open(snakemake.input['matrix_uk'], 'rb') as fr:
    matrix_uk = pickle.load(fr)

with open(snakemake.input['matrix_bo'], 'rb') as fr:
    matrix_bo = pickle.load(fr)

with open(snakemake.input['file_data_uk'], 'rb') as fr:
    data_uk = pickle.load(fr)

with open(snakemake.input['file_data_bo'], 'rb') as fr:
    data_bo = pickle.load(fr)
    

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

#=============================================
# First figure

fig, axs = plt.subplots(
    2, 2, figsize=(12, 8),
    height_ratios=[1, 0.7],
)
fig.tight_layout(h_pad=2, rect=[0.02, 0.04, 1, 0.98])

tri_lo = np.tril_indices(matrix_uk['straps_cov_nc'][1].shape[0], k=0)
for i in [0, 1, 2]: 
    matrix_uk['straps_cov_nc'][i][tri_lo] = np.nan

tri_lo = np.tril_indices(matrix_bo['straps_cov_nc'][1].shape[0], k=0)
for i in [0, 1, 2]: 
	matrix_bo['straps_cov_nc'][i][tri_lo] = np.nan

k, l = (0, 0)
labels = [f"$\\Delta p_{i}$" for i, t in enumerate(matrix_uk['times'][:-1])]
ac.plot_covmat_ci(
    matrix_uk['straps_cov_nc'],
    axs[k, l],
    delta_labels=labels,
    cbar_kws={
        'label': 'covariance',
        'format': formatter,
        'location': 'left',
    },
)
axs[k, l].set_title("A", loc='left', fontdict={'fontweight': 'bold'})
axs[k, l].set_title('UK', fontweight='bold')
axs[k, l].set(xticklabels=[])

k, l = (0, 1)
labels = [f"$\\Delta p_{i}$" for i, t in enumerate(matrix_bo['times'][:-1])]
ac.plot_covmat_ci(
    matrix_bo['straps_cov_nc'],
    axs[k, l],
    delta_labels=labels,
    cbar_kws={
        'label': 'covariance',
        'format': formatter,
        'location': 'left',
    },
)
axs[k, l].set_title("B", loc='left', fontdict={'fontweight': 'bold'})
axs[k, l].set_title('Bohemia', fontweight='bold')
axs[k, l].set(xticklabels=[])

fmts = ['-o', '-s', '-^']
labels = ['WHG-like', 'EEF-like', 'Steppe-like']
delta_list_uk = [f"$\\Delta p_{{{int(t)}}}$" for t in range(len(data_uk['times']) - 1)]
k, l = (1, 0)
qorder = [0,1,2]
for i in range(len(labels)):
    axs[k, l].plot(data_uk['times'], data_uk['Q'][:,qorder[i]], fmts[i], label=labels[i], color=colors_oi[i])
    # Q order WHG, EEF, Steppe
for x1, x2, txt in zip(data_uk['times'][:-1], data_uk['times'][1:], delta_list_uk):
    _ = axs[k, l].text(x2+(x1 - x2)/2, 0.95, txt, ha='center')
for i, t in enumerate(data_uk['times']):
    _ = axs[k, l].text(t, 0.85, str(i), ha='center')
for x1, x2 in zip(data_uk['times'][1::2], data_uk['times'][2::2]):
    _ = axs[k, l].axvspan(x1, x2, facecolor='grey', alpha=0.10)
axs[k, l].set_xlim(data_uk['times'][0] + time_padding, data_uk['times'][-1] - time_padding)
axs[k, l].set_ylim(top=1)
axs[k, l].set_ylabel("Mean ancestry proportion")
axs[k, l].set_xlabel("Time (years BP)")
# axs[k, l].legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=3)
axs[k, l].set_title("C", loc='left', fontdict={'fontweight': 'bold'})

delta_list_bo = [f"$\\Delta p_{{{int(t)}}}$" for t in range(len(data_bo['times']) - 1)]
k, l = (1, 1)
qorder = [1,0,2]
for i in range(len(labels)):
    axs[k, l].plot(data_bo['times'], data_bo['Q'][:,qorder[i]], fmts[i], label=labels[i], color=colors_oi[i])
    # Q order EEF, WHG, Steppe
for x1, x2, txt in zip(data_bo['times'][:-1], data_bo['times'][1:], delta_list_bo):
    _ = axs[k, l].text(x2+(x1 - x2)/2, 0.95, txt, ha='center')
for i, t in enumerate(data_bo['times']):
    _ = axs[k, l].text(t, 0.85, str(i), ha='center')
for x1, x2 in zip(data_bo['times'][1::2], data_bo['times'][2::2]):
    _ = axs[k, l].axvspan(x1, x2, facecolor='grey', alpha=0.10)
axs[k, l].set_xlim(data_bo['times'][0] + time_padding, data_bo['times'][-1] - time_padding)
axs[k, l].set_ylim(top=1)
axs[k, l].set_ylabel("Mean ancestry proportion")
axs[k, l].set_xlabel("Time (years BP)")
# axs[k, l].legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=3)
axs[k, l].set_title("D", loc='left', fontdict={'fontweight': 'bold'})

handles, labels = axs[1, 1].get_legend_handles_labels()
fig.legend(
    handles, labels,
    loc='lower center', bbox_to_anchor=(0.5, 0), ncol=3
)

# add connections
from matplotlib.patches import ConnectionPatch

for i, t in enumerate(data_uk['times']):
    con = ConnectionPatch(
        xyA=(i, 6), coordsA=axs[0,0].transData,
        xyB=(t, 1), coordsB=axs[1,0].transData,
        alpha=0.6,
    )
    fig.add_artist(con)
for i, t in enumerate(data_bo['times']):
    con = ConnectionPatch(
        xyA=(i, 6), coordsA=axs[0,1].transData,
        xyB=(t, 1), coordsB=axs[1,1].transData,
        alpha=0.6,
    )
    fig.add_artist(con)

fig.savefig(snakemake.output['fig_data_matrix_ancestry'])


#=============================================
# Second figure

fig, axs = plt.subplots(
    4, 2, figsize=(12, 12), layout='tight',
    height_ratios=[0.08, 1, 1, 1],
)

k, l = (0, 0)
for i, t1, t2, txt in zip(range(7), data_uk['times'][:-1], data_uk['times'][1:], delta_list_uk):
    _ = axs[k, l].axvspan(t1, t2, facecolor=colors_oi[i])
    _ = axs[k, l].text(t2+(t1 - t2)/2, 0.5, txt, ha='center', va='center', color='white', fontdict={'fontweight': 'bold'})
for i, t in enumerate(data_uk['times']):
    _ = axs[k, l].text(t, 1.1, str(i), ha='center')
axs[k, l].set(yticklabels=[])
axs[k, l].tick_params(left=False)
axs[k, l].set_ylabel("$\\Delta p_i$", rotation=0, verticalalignment='center')
axs[k, l].set_xlim(data_uk['times'][0], data_uk['times'][-1])
axs[k, l].set_xlabel("Time (years BP)")
axs[k, l].set_title("UK", fontweight='bold', pad=15)
axs[k, l].set_xticks(
    ticks=data_uk['times'], labels=[int(t) for t in data_uk['times']],
    rotation=45, horizontalalignment='right', rotation_mode='anchor',
)
# axs[k, l].tick_params(pad=0)

k, l = (0, 1)
for i, t1, t2, txt in zip(range(7), data_bo['times'][:-1], data_bo['times'][1:], delta_list_bo):
    _ = axs[k, l].axvspan(t1, t2, facecolor=colors_oi[i])
    _ = axs[k, l].text(t2+(t1 - t2)/2, 0.5, txt, ha='center', va='center', color='white', fontdict={'fontweight': 'bold'})
for i, t in enumerate(data_bo['times']):
    _ = axs[k, l].text(t, 1.1, str(i), ha='center')
axs[k, l].set(yticklabels=[])
axs[k, l].tick_params(left=False)
axs[k, l].set_ylabel("$\\Delta p_i$", rotation=0, verticalalignment='center')
axs[k, l].set_xlim(data_bo['times'][0], data_bo['times'][-1] - 50)
axs[k, l].set_xlabel("Time (years BP)")
axs[k, l].set_title('Bohemia', fontweight='bold', pad=15)
axs[k, l].set_xticks(
    ticks=data_bo['times'], labels=[int(t) for t in data_bo['times']],
    rotation=45, horizontalalignment='right', rotation_mode='anchor',
)
# axs[k, l].tick_params(pad=0)
axs[k, l].set_facecolor('grey')

# uk

x_shift = 0.1
new_times = np.array(range(len(data_uk['times'])))
k, l = (1, 0)
ac.cov_lineplot(new_times, data_uk['straps_cov_nc'], axs[k, l], colors=colors_oi, d=x_shift, labels=delta_list_uk)
axs[k, l].set_xlim(new_times[1] - x_shift, new_times[-2] + 3 * x_shift)
axs[k, l].hlines(y=0, xmin=0, xmax=new_times[-1] + 3 * x_shift, linestyles='dotted', colors='grey')
axs[k, l].set_ylabel("Cov($\\Delta p_i$, $\\Delta p_j$)")
axs[k, l].set_xlabel("$\\Delta p_j$")
axs[k, l].set_title('Before admixture correction')
axs[k, l].set_title("A", loc='left', fontdict={'fontweight': 'bold'})
axs[k, l].xaxis.set_major_locator(loc)
axs[k, l].yaxis.set_major_formatter(formatter)
axs[k, l].set_xticks(ticks=range(1, 6), labels=delta_list_uk[1:])

k, l = (2, 0)
ac.cov_lineplot(new_times, data_uk['straps_cov'], axs[k, l], colors=colors_oi, d=x_shift, labels=delta_list_uk, ylim=axs[1, 0].get_ylim())
axs[k, l].set_xlim(new_times[1] - x_shift, new_times[-2] + 3 * x_shift)
axs[k, l].hlines(y=0, xmin=0, xmax=new_times[-1] + 3 * x_shift, linestyles='dotted', colors='grey')
axs[k, l].set_ylabel("Cov($\\Delta p_i$, $\\Delta p_j$)")
axs[k, l].set_xlabel('$\\Delta p_j$')
axs[k, l].set_title('After admixture correction')
axs[k, l].set_title("C", loc='left', fontdict={'fontweight': 'bold'})
axs[k, l].xaxis.set_major_locator(loc)
axs[k, l].yaxis.set_major_formatter(formatter)
axs[k, l].set_xticks(ticks=range(1, 6), labels=delta_list_uk[1:])

k, l = (3, 0)
ac.plot_ci_line(new_times[1:] + x_shift, np.stack(data_uk['straps_G_nc']).T, ax=axs[k, l], linestyle='dashed', marker='o', label='$G_{nc}$')
ac.plot_ci_line(new_times[1:], np.stack(data_uk['straps_G']).T, ax=axs[k, l], marker='o', label='$G$')
ac.plot_ci_line(new_times[1:] - x_shift, np.stack(data_uk['straps_Ap']).T, ax=axs[k, l], color='blue', marker='s', label='$A$')
axs[k, l].set_xlim(new_times[1] - 2*x_shift, new_times[-1] + 2*x_shift)
axs[k, l].hlines(y=0, xmin=new_times[-1], xmax=new_times[1], colors='grey', linestyles='dotted')
axs[k, l].set_ylim(ymax=1)
axs[k, l].set_xlabel('$p_t$')
axs[k, l].set_ylabel("Proportion of variance ($p_t - p_{0}$)")
axs[k, l].set_title("E", loc='left', fontdict={'fontweight': 'bold'})
axs[k, l].xaxis.set_major_locator(loc)
axs[k, l].set_xticks(ticks=range(1, 7), labels=[f'$p_{i}$' for i in range(1, 7)])
axs[k, l].set_title('Variance decomposition')
for ci, t in zip(data_uk['straps_G'], new_times[1:]):
    if ci[0]*ci[2] > 0:
        axs[k, l].annotate("*", xy=(t, 0.2))
axs[k, l].legend(loc='upper right')


# bo

x_shift = 0.1
new_times = np.array(range(len(data_bo['times'])))
k, l = (1, 1)
ac.cov_lineplot(new_times, data_bo['straps_cov_nc'], axs[k, l], colors=colors_oi, d=x_shift, labels=delta_list_bo)
axs[k, l].set_xlim(new_times[1] - x_shift, new_times[-2] + 3 * x_shift)
axs[k, l].hlines(y=0, xmin=0, xmax=new_times[-1] + 3 * x_shift, linestyles='dotted', colors='grey')
axs[k, l].set_ylabel("Cov($\\Delta p_i$, $\\Delta p_j$)")
axs[k, l].set_xlabel('$\\Delta p_j$')
axs[k, l].set_title('Before admixture correction')
axs[k, l].set_title("B", loc='left', fontdict={'fontweight': 'bold'})
axs[k, l].xaxis.set_major_locator(loc)
axs[k, l].yaxis.set_major_formatter(formatter)
axs[k, l].set_xticks(ticks=range(1, 6), labels=delta_list_bo[1:])

k, l = (2, 1)
ac.cov_lineplot(new_times, data_bo['straps_cov'], axs[k, l], colors=colors_oi, d=x_shift, labels=delta_list_bo, ylim=axs[1, 1].get_ylim())
axs[k, l].set_xlim(new_times[1] - x_shift, new_times[-2] + 3 * x_shift)
axs[k, l].hlines(y=0, xmin=0, xmax=new_times[-1] + 3 * x_shift, linestyles='dotted', colors='grey')
axs[k, l].set_ylabel("Cov($\\Delta p_i$, $\\Delta p_j$)")
axs[k, l].set_xlabel('$\\Delta p_j$')
axs[k, l].set_title('After admixture correction')
axs[k, l].set_title("D", loc='left', fontdict={'fontweight': 'bold'})
axs[k, l].xaxis.set_major_locator(loc)
axs[k, l].yaxis.set_major_formatter(formatter)
axs[k, l].set_xticks(ticks=range(1, 6), labels=delta_list_bo[1:])

k, l = (3, 1)
ac.plot_ci_line(new_times[1:] + x_shift, np.stack(data_bo['straps_G_nc']).T, ax=axs[k, l], linestyle='dashed', marker='o', label='$G_{nc}$')
ac.plot_ci_line(new_times[1:], np.stack(data_bo['straps_G']).T, ax=axs[k, l], marker='o', label='$G$')
ac.plot_ci_line(new_times[1:] - x_shift, np.stack(data_bo['straps_Ap']).T, ax=axs[k, l], color='blue', marker='s', label='$A$')
axs[k, l].set_xlim(new_times[1] - 2*x_shift, new_times[-1] + 2*x_shift)
axs[k, l].hlines(y=0, xmin=new_times[-1], xmax=new_times[1], colors='grey', linestyles='dotted')
axs[k, l].set_ylim(ymax=1)
axs[k, l].set_xlabel('$p_t$')
axs[k, l].set_ylabel("Proportion of variance ($p_t - p_{0}$)")
axs[k, l].set_title("F", loc='left', fontdict={'fontweight': 'bold'})
axs[k, l].xaxis.set_major_locator(loc)
axs[k, l].set_xticks(ticks=range(1, 7), labels=[f'$p_{i}$' for i in range(1, 7)])
axs[k, l].set_title('Variance decomposition')
for ci, t in zip(data_bo['straps_G'], new_times[1:]):
    if ci[0]*ci[2] > 0:
        axs[k, l].annotate("*", xy=(t, 0.2))


# handles, labels = axs[3, 0].get_legend_handles_labels()
# fig.legend(
#     handles, labels,
#     loc='lower center', bbox_to_anchor=(0.5, 0.27), ncol=3
# )

fig.savefig(snakemake.output['fig_data_covlines'])