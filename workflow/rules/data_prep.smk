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
		multiext("data/v50.0_1240k_public", '.geno', '.snp', '.ind')
	shell:
		"""
		cd data
		tar -xf v50.0_1240K_public.tar
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
		multiext("data/v50.0_1240k_public", '.geno', '.snp', '.ind'),
		multiext("data/brit", '.geno', '.snp', '.ind'),
	output:
		"data/mergit_Patterson2022.par",
		"data/v50.0_1240k_public_filtered.ind",
	conda: "../envs/analyses.yaml"
	script:
		"../scripts/data_filter_merge.R"

rule data_filter_merge_2:
	input:
		"data/mergit_Patterson2022.par",
		multiext("data/v50.0_1240k_public", '.geno', '.snp', '.ind'),
		multiext("data/brit", '.geno', '.snp', '.ind'),
	output:
		multiext("data/Patterson2022", ".bed", ".bim", ".fam"),
	params:
		prefix = lambda w, output: output[0].replace(".bed", ""),
		tmp = lambda w, output: output[0].replace("data/", "data/tmp_").replace(".bed", ""),
	conda: "../envs/analyses.yaml"
	shell:
		"""
		mergeit -p {input[0]}
		plink --bfile {params.tmp} \
		--make-bed --allow-no-sex --out {params.prefix} && \
		rm {params.tmp}.*
		"""

rule convert_plink2sgkit:
	input:
		multiext("data/Patterson2022", ".bed", ".bim", ".fam"),
		metadata = "data/Patterson2022_TableS3.tsv",
	output:
		zarr_store = directory("data/Patterson2022.zarr"),
		meta_in = "data/Patterson2022_meta_in.tsv",
		meta_out = "data/Patterson2022_meta_out.tsv",
	params:
		path = lambda w, input: input[0].replace(".bed", ""),
	conda: "../envs/analyses.yaml"
	script:
		"../scripts/convert_plink2sgkit.py"
