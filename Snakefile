# Snakemake workflow

configfile: "config.yaml"

include: "workflow/rules/data_prep_Patterson2022.smk"
include: "workflow/rules/data_prep_Papac2021.smk"
include: "workflow/rules/smartpca.smk"

rule all:
	input:
		rules.convert_plink2sgkit_Patterson2022.output,
		rules.smartpca_Patterson2022.output,
		rules.smartpca_Patterson2022_proj.output,
		rules.prepare_maps.output,
		rules.make_bed_Papac2021.output,
