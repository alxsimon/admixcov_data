# Snakemake workflow

configfile: "config.yaml"

include: "workflow/rules/data_prep.smk"

rule all:
	input:
		rules.data_filter_merge_2.output