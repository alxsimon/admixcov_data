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
		multiext("data/Patterson2022/brit", '.geno', '.snp', '.ind')
	params: 
		url = config['Patterson2022_url']
	shell:
		"""
		wget {params.url} -O "data/Patterson2022/Patterson2022.zip"
		cd data/Patterson2022/
		unzip Patterson2022.zip && rm Patterson2022.zip
		rm -r __MACOSX
		"""

rule download_Fowler2021:
	output:
		multiext("data/Patterson2022/Hazleton_v51", '.geno', '.snp', '.ind')
	params: 
		url = config['Fowler2021_url']
	shell:
		"""
		wget {params.url} -O "data/Patterson2022/Fowler2021.zip"
		cd data/Patterson2022/
		unzip Fowler2021.zip && rm Fowler2021.zip
		"""

rule data_filter_merge_1:
	input: 
		sample_list = "data/Patterson2022/list_selected_samples.txt",
		aadr_ind = "data/v50.0_1240k_public.ind",
		brit_ind = "data/Patterson2022/brit.ind",
		hazleton_ind = "data/Patterson2022/Hazleton_v51.ind",
		aadr_modern = "data/list_samples_aadr_modern.txt",
	output:
		"data/v50.0_1240k_public_filtered_Patterson2022.ind",
		"data/Patterson2022/brit_filtered.ind",
		"data/Patterson2022/Hazleton_v51_filtered.ind",
	conda: "../envs/r-env.yaml"
	script:
		"../scripts/data_filter_merge_Patterson2022.R"


rule merge_Patterson2022_1:
	input:
		multiext("data/Patterson2022/brit", ".geno", ".snp", "_filtered.ind"),
		multiext("data/Patterson2022/Hazleton_v51", ".geno", ".snp", "_filtered.ind"),
	output:
		"data/Patterson2022/mergit_1.par",
		temp(multiext("data/Patterson2022/tmp_1", ".geno", ".snp", ".ind")),
	log:
		"logs/mergit_Patterson2022_1.log"
	conda: "../envs/eigensoft.yaml"
	shell:
		"""
		cat << EOF > {output[0]}
geno1: {input[0]}
snp1: {input[1]}
ind1: {input[2]}
geno2: {input[3]}
snp2: {input[4]}
ind2: {input[5]}
outputformat: EIGENSTRAT
genooutfilename: {output[1]}
snpoutfilename: {output[2]}
indoutfilename: {output[3]}
allowdups: NO
hashcheck: NO
EOF

		mergeit -p {output[0]} > {log} 2>&1
		"""

rule merge_Patterson2022_2:
	input:
		multiext("data/Patterson2022/tmp_1", ".geno", ".snp", ".ind"),
		multiext("data/v50.0_1240k_public", ".geno", ".snp", "_filtered_Patterson2022.ind"),
	output:
		"data/Patterson2022/mergit_2.par",
		temp(multiext("data/Patterson2022/tmp_2", ".bed", ".bim", ".fam")),
	log:
		"logs/mergit_Patterson2022_2.log"
	conda: "../envs/eigensoft.yaml"
	shell:
		"""
		cat << EOF > {output[0]}
geno1: {input[0]}
snp1: {input[1]}
ind1: {input[2]}
geno2: {input[3]}
snp2: {input[4]}
ind2: {input[5]}
outputformat: PACKEDPED
genooutfilename: {output[1]}
snpoutfilename: {output[2]}
indoutfilename: {output[3]}
allowdups: NO
hashcheck: NO
EOF

		mergeit -p {output[0]} > {log} 2>&1
		"""


rule make_bed_Patterson2022:
	input:
		multiext("data/Patterson2022/tmp_2", ".bed", ".bim", ".fam"),
	output:
		multiext("data/Patterson2022/Patterson2022", ".bed", ".bim", ".fam"),
	params:
		prefix_in = lambda w, input: input[0].replace(".bed", ""),
		prefix_out = lambda w, output: output[0].replace(".bed", ""),
	conda: "../envs/eigensoft.yaml"
	shell:
		"""
		plink --bfile {params.prefix_in} \
		--make-bed --allow-no-sex --out {params.prefix_out}
		"""

rule prepare_maps:
	input:
		bmap_files = expand('data/Murphy2021_Bvalues/CADD_bestfit/chr{num}.bmap.txt', num=range(1, 23)),
		rmap_files = expand('data/Bherer2017_Refined_EUR_genetic_map_b37/sexavg_chr{num}.txt', num=range(1, 23)),
		chr_len = 'data/hg19_chr_len.txt',
	output:
		bmap_out = 'data/Murphy2021_Bvalues_compiled.bmap.txt',
		rmap_out = 'data/Bherer2017_Refined_EUR_genetic_map_sexavg.rmap.txt',
	conda: "../envs/py-env.yaml"
	script:
		"../scripts/prepare_maps.py"

rule convert_plink2sgkit_Patterson2022:
	input:
		multiext("data/Patterson2022/Patterson2022", ".bed", ".bim", ".fam"),
		metadata_S3 = "data/Patterson2022/Patterson2022_TableS3.tsv",
		metadata_S5 = "data/Patterson2022/Patterson2022_TableS5_corrected.tsv",
		modern_anno = "data/v50.0_1240k_public_modern_selection.anno",
		bval_table = 'data/Murphy2021_Bvalues_compiled.bmap.txt',
		rval_table = 'data/Bherer2017_Refined_EUR_genetic_map_sexavg.rmap.txt',
	output:
		zarr_store = directory("data/Patterson2022/Patterson2022.zarr"),
		meta_out = "data/Patterson2022/Patterson2022_metadata.tsv",
		out_ind = "data/Patterson2022/Patterson2022.ind",
	params:
		path = lambda w, input: input[0].replace(".bed", ""),
	conda: "../envs/py-env.yaml"
	script:
		"../scripts/convert_plink2sgkit_Patterson2022.py"
