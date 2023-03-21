# Snakemake workflow

configfile: "config.yaml"

include: "workflow/rules/data_prep_Patterson2022.smk"
include: "workflow/rules/data_prep_Papac2021.smk"
include: "workflow/rules/smartpca.smk"
include: "workflow/rules/analyses.smk"

rule all:
	input:
		rules.smartpca_Patterson2022.output,
		rules.smartpca_Patterson2022_proj.output,
		rules.prepare_maps.output,
		'results/Patterson2022/main_figure_Patterson2022.pdf',
		'results/Papac2021/main_figure_Papac2021.pdf',
