#!python3

import pickle
import matplotlib.pyplot as plt

import matplotlib.ticker as tkr
loc = tkr.MultipleLocator(base=1.0)
# sci notation formatter
formatter = tkr.ScalarFormatter(useMathText=True)
formatter.set_scientific(True)
formatter.set_powerlimits((0, 0))

with open(snakemake.input['file_data_uk_ascert'], 'rb') as fr:
    data_uk_ascert = pickle.load(fr)
with open(snakemake.input['file_data_bo_ascert'], 'rb') as fr:
    data_bo_ascert = pickle.load(fr)

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

fig, axs = plt.subplots(1, 2, figsize=(12, 4), layout="constrained")

y_cat = [
	"13-Denisova-San",
	"5-Yoruba",
	"4-San",
	"3-Papuan1",
	"2-Han",
    "1-French",
    "total",
]
# UK subplot
data = data_uk
data_ascert = data_uk_ascert
for val, fmt, col, label in zip(
    ['straps_G', 'straps_Ap'],
    ['o', 's'],
    ['k', 'b'],
    ["$G$", "$A$"],
):
    x_err = [[], []]
    x_mean = []
    y_N_snp = []

    for i in list(range(6))[-1::-1]:
        x_mean.append(
            data_ascert[i][val][-1][1]
        )
        x_err[0].append(
            x_mean[-1] - data_ascert[i][val][-1][0]
        )
        x_err[1].append(
            data_ascert[i][val][-1][2] - x_mean[-1]
        )
        y_N_snp.append(data_ascert[i]["N_snp"])
    x_mean.append(
        data[val][-1][1]
    )
    x_err[0].append(
        x_mean[-1] - data[val][-1][0]
    )
    x_err[1].append(
        data[val][-1][2] - x_mean[-1]
    )
    y_N_snp.append(data["N_snp"])

    axs[0].errorbar(
        x=x_mean,
        y=[cat + "\n" + str(N) for cat, N in zip(y_cat, y_N_snp)],
        xerr=x_err,
        fmt=fmt,
        color=col,
        label=label,
    )


# Bohemia subplot
data = data_bo
data_ascert = data_bo_ascert
for val, fmt, col, label in zip(
    ['straps_G', 'straps_Ap'],
    ['o', 's'],
    ['k', 'b'],
    ["$G$", "$A$"],
):
    x_err = [[], []]
    x_mean = []
    y_N_snp = []

    for i in list(range(6))[-1::-1]:
        x_mean.append(
            data_ascert[i][val][-1][1]
        )
        x_err[0].append(
            x_mean[-1] - data_ascert[i][val][-1][0]
        )
        x_err[1].append(
            data_ascert[i][val][-1][2] - x_mean[-1]
        )
        y_N_snp.append(data_ascert[i]["N_snp"])
    x_mean.append(
        data[val][-1][1]
    )
    x_err[0].append(
        x_mean[-1] - data[val][-1][0]
    )
    x_err[1].append(
        data[val][-1][2] - x_mean[-1]
    )
    y_N_snp.append(data["N_snp"])

    axs[1].errorbar(
        x=x_mean,
        y=[cat + "\n" + str(N) for cat, N in zip(y_cat, y_N_snp)],
        xerr=x_err,
        fmt=fmt,
        color=col,
        label=label,
    )

axs[0].set_title("A", loc='left', fontdict={'fontweight': 'bold'})
axs[0].set_title('UK', fontweight='bold')
axs[0].set_xlabel("Proportion of variance")
axs[1].set_title("B", loc='left', fontdict={'fontweight': 'bold'})
axs[1].set_title('Bohemia', fontweight='bold')
axs[1].set_xlabel("Proportion of variance")

axs[0].vlines(x=0, ymin=0, ymax=6, linestyles='dotted', colors='grey')
axs[1].vlines(x=0, ymin=0, ymax=6, linestyles='dotted', colors='grey')

fig.savefig(snakemake.output['fig'])
