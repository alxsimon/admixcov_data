# Snakemake workflow

configfile: "config.yaml"

include: "workflow/rules/data_prep.smk"

rule all:
	input:
		rules.convert_plink2sgkit.output