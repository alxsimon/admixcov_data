rule download_AADR:
	output:
		multiext("data/v50.0_1240K_public/v50.0_1240k_public", '.geno', '.snp', '.ind')
	params:
		url = config['AADR_url']
	shell:
		"""
		wget {params.url} -O "data/v50.0_1240K_public.tar"
		cd data
		tar -x v50.0_1240K_public.tar && rm v50.0_1240K_public.tar
		"""

rule download_Patterson2022:
	output:
		multiext("data/brit", '.geno', '.snp', '.ind')
	params: 
		url = config['Patterson2022_url']
	shell:
		"""
		wget {params.url} -O "data/tmp_data.zip"
		cd data
		unzip tmp_data.zip && rm tmp_data.zip
		"""

rule data_filter_merge_1:
	input: 
		"data/Patterson2022_TableS3.tsv",
		multiext("data/v50.0_1240K_public/v50.0_1240k_public", '.geno', '.snp', '.ind'),
		multiext("data/brit", '.geno', '.snp', '.ind'),
	output:
		"data/mergit_Patterson2022.par",
		"data/v50.0_1240K_public/v50.0_1240k_public_filtered.ind",
	conda: "../envs/analyses.yaml"
	script:
		"scripts/data_filter_merge.R"

rule data_filter_merge_2:
	input:
		"data/mergit_Patterson2022.par",
		multiext("data/v50.0_1240K_public/v50.0_1240k_public", '.geno', '.snp', '.ind'),
		multiext("data/brit", '.geno', '.snp', '.ind'),
	output:
		multiext("data/Patterson2022", ".ped", ".map", ".fam"),
	conda: "../envs/analyses.yaml"
	shell:
		"""
		mergeit -p {input[0]}
		"""