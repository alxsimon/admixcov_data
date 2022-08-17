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

rule download_Fowler2021:
	output:
		multiext("data/Hazleton_v51", '.geno', '.snp', '.ind')
	params: 
		url = config['Fowler2021_url']
	shell:
		"""
		wget {params.url} -O "data/tmp_data.zip"
		cd data
		unzip tmp_data.zip && rm tmp_data.zip
		"""

rule data_filter_merge_1:
	input: 
		"data/list_selected_samples.txt",
		"data/v50.0_1240k_public.ind",
		"data/brit.ind",
		"data/Hazleton_v51.ind",
	output:
		"data/mergit_intermediate.par",
		"data/mergit_Patterson2022.par",
		"data/v50.0_1240k_public_filtered.ind",
		"data/brit_filtered.ind",
		"data/Hazleton_v51_filtered.ind",
	conda: "../envs/analyses.yaml"
	script:
		"../scripts/data_filter_merge.R"

rule data_filter_merge_2:
	input:
		"data/mergit_intermediate.par",
		"data/mergit_Patterson2022.par",
		multiext("data/v50.0_1240k_public", '.geno', '.snp', '_filtered.ind'),
		multiext("data/brit", '.geno', '.snp', '_filtered.ind'),
		multiext("data/Hazleton_v51", '.geno', '.snp', '_filtered.ind'),
	output:
		multiext("data/Patterson2022", ".bed", ".bim", ".fam"),
	params:
		prefix = lambda w, output: output[0].replace(".bed", ""),
	conda: "../envs/analyses.yaml"
	shell:
		"""
		mergeit -p {input[0]}
		mergeit -p {input[1]}
		plink --bfile data/tmp_2 \
		--make-bed --allow-no-sex --out {params.prefix} && \
		rm data/tmp_*
		"""

rule prepare_maps:
	input:
		bmap_files = expand('data/Murphy2021_Bvalues/CADD_bestfit/chr{num}.bmap.txt', num=range(1, 23)),
		rmap_files = expand('data/Bherer2017_Refined_EUR_genetic_map_b37/sexavg_chr{num}.txt', num=range(1, 23)),
		chr_len = 'data/hg19_chr_len.txt',
	output:
		bmap_out = 'data/Murphy2021_Bvalues_compiled.bmap.txt',
		rmap_out = 'data/Bherer2017_Refined_EUR_genetic_map_sexavg.rmap.txt',
	conda: "../envs/analyses.yaml"
	script:
		"../scripts/prepare_maps.py"

rule convert_plink2sgkit:
	input:
		multiext("data/Patterson2022", ".bed", ".bim", ".fam"),
		metadata_S3 = "data/Patterson2022_TableS3.tsv",
		metadata_S5 = "data/Patterson2022_TableS5_corrected.tsv",
		bval_table = 'data/Murphy2021_Bvalues_compiled.bmap.txt',
		rval_table = 'data/Bherer2017_Refined_EUR_genetic_map_sexavg.rmap.txt',
	output:
		zarr_store = directory("data/Patterson2022.zarr"),
		meta_out = "data/Patterson2022_metadata.tsv",
		out_ind = "data/Patterson2022.ind",
	params:
		path = lambda w, input: input[0].replace(".bed", ""),
	conda: "../envs/analyses.yaml"
	script:
		"../scripts/convert_plink2sgkit.py"
