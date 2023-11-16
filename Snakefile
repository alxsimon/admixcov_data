# Snakemake workflow

configfile: "config.yaml"

include: "workflow/rules/data_prep_Patterson2022.smk"
include: "workflow/rules/data_prep_Papac2021.smk"
include: "workflow/rules/smartpca.smk"
include: "workflow/rules/analyses.smk"
include: "workflow/rules/ascertainment.smk"

rule all:
	input:
		"results.tar.gz",
		rules.analysis_papac.output,
		rules.main_figures.output,
		expand(
			"results/{dataset}/ascertainment_info/Rohland2022_supp_chr{chr}.txt",
			dataset=["Patterson2022", "Papac2021"],
			chr=range(1, 22),
		)


rule download_AADR:
	output:
		temp("data/v50.0_1240K_public.tar")
	params:
		url = config['AADR_url']
	shell:
		"""
		wget {params.url} -O {output}
		"""

rule extract_AADR:
	input:
		"data/v50.0_1240K_public.tar"
	output:
		multiext("data/v50.0_1240k_public", '.geno', '.snp', '.ind', '.anno')
	shell:
		"""
		cd data
		tar -xf v50.0_1240K_public.tar
		"""


rule prepare_maps:
	input:
		bmap_files = expand('data/Murphy2021_Bvalues/CADD_bestfit/chr{num}.bmap.txt', num=range(1, 23)),
		rmap_files = expand('data/Bherer2017_Refined_EUR_genetic_map_b37/sexavg_chr{num}.txt', num=range(1, 23)),
		chr_len = 'data/hg19_chr_len.txt',
	output:
		bmap_out = 'data/Murphy2021_Bvalues_compiled.bmap.txt',
		rmap_out = 'data/Bherer2017_Refined_EUR_genetic_map_sexavg.rmap.txt',
	conda: "workflow/envs/py-env.yaml"
	script:
		"workflow/scripts/prepare_maps.py"


rule archive:
	input:
		rules.analysis_patterson.output,
		rules.analysis_papac.output,
		rules.smartpca_Patterson2022.output,
		rules.smartpca_Patterson2022_proj.output,
		rules.prepare_maps.output,
		rules.analysis_patterson.output,
		rules.analysis_papac.output,
		rules.main_figures.output,
		rules.analysis_patterson_split.output,
	output:
		"results.tar.gz"
	shell:
		"tar -czf {output} results"