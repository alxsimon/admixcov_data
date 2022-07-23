# Snakemake workflow

configfile: "config.yaml"

include: "workflow/rules/data_prep.smk"
include: "workflow/rules/smartpca.smk"

rule all:
	input:
		rules.convert_plink2sgkit.output,
		rules.smartpca.output,
		rules.smartpca_ref.output,