# Snakemake workflow

configfile: "config.yaml"

include: "workflow/rules/data_prep_Patterson2022.smk"
include: "workflow/rules/data_prep_Papac2021.smk"
include: "workflow/rules/smartpca.smk"
include: "workflow/rules/analyses.smk"

rule archive:
	input:
		rules.analysis_patterson.output,
		rules.analysis_papac.output,
		rules.smartpca_Patterson2022.output,
		rules.smartpca_Patterson2022_proj.output,
		rules.prepare_maps.output,
		rules.analysis_patterson.output,
		rules.analysis_papac.output,
		rules.figure_matrices.output,
		rules.analysis_patterson_split.output,
	output:
		"results.tar.gz"
	shell:
		"tar -czf {output} results"

rule all:
	input:
		rules.archive.output,