rule grab_common_snps:
	input:
		used_snps = "results/{dataset}/used_snps_{dataset}.tsv",
		rohland_chr = "data/Rohland2022/supp_gr.276728.122_Supplemental_Data_7_chr{chr}.txt.gz",
	output:
		out = "results/{dataset}/ascertainment_info/Rohland2022_supp_chr{chr}.txt",
	conda:
		"../envs/command-line.yaml"
	resources:
		mem = "3G"
	shell:
		"""
		awk 'NR==FNR {{a[$1]; b[$2]; next}} ($1 in a) && ($2 in b)' \
		<(awk -v chr="{wildcards.chr}" '$1==chr' {input.used_snps}) \
		<(zcat {input.rohland_chr}) \
		> {output.out}
		"""


rule Patterson2022_ascert:
	input:
		zarr = 'data/Patterson2022/Patterson2022.zarr',
		ascert_info = expand(
			"results/Patterson2022/ascertainment_info/Rohland2022_supp_chr{chr}.txt",
			chr=range(1, 22),
		),
	output:
		fig_data = 'results/Patterson2022/fig_data_Patterson2022_ascert.pickle',
	conda:
		"../envs/py-env.yaml"
	script:
		'../scripts/analysis_patterson_ascert.py'


rule Papac2021_ascert:
	input:
		zarr = 'data/Papac2021/Papac2021.zarr',
		ascert_info = expand(
			"results/Papac2021/ascertainment_info/Rohland2022_supp_chr{chr}.txt",
			chr=range(1, 22),
		),
	output:
		fig_data = 'results/Papac2021/fig_data_Papac2021_ascert.pickle',
	conda:
		"../envs/py-env.yaml"
	script:
		'../scripts/analysis_papac_ascert.py'