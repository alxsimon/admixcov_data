
rule download_Papac2021:
	output:
		multiext("data/Papac2021/Papac2021_forUpload", '.geno', '.snp', '.ind')
	params: 
		url = config['Papac2021_url']
	shell:
		"""
		wget {params.url} -O "data/Papac2021/Papac2021.zip"
		cd data/Papac2021/
		unzip Papac2021.zip && rm Papac2021.zip
		"""

rule data_filter_merge_Papac2021:
	input: 
		multiext("data/v50.0_1240k_public", ".geno", ".snp", ".ind"),
		multiext("data/Papac2021/Papac2021_forUpload", '.geno', '.snp', '.ind'),
		"data/Papac2021/Modern_samples_Papac2021.ind",
		"data/Papac2021/Papac2021_TableS5_ancient_published.tsv",
	output:
		"data/Papac2021/Papac2021_forUpload.modif.snp",
		"data/v50.0_1240k_public_filtered_Papac2021.ind",
	conda: "../envs/r-env.yaml"
	script:
		"../scripts/data_filter_merge_Papac2021.R"


rule merge_Papac2021:
	input:
		multiext("data/Papac2021/Papac2021_forUpload", ".geno", ".modif.snp", ".ind"),
		multiext("data/v50.0_1240k_public", ".geno", ".snp", "_filtered_Papac2021.ind"),
	output:
		"data/Papac2021/mergit.par",
		temp(multiext("data/Papac2021/tmp", ".bed", ".bim", ".fam")),
	log:
		"logs/mergit_Papac2021.log"
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
allowdups: YES
hashcheck: NO
EOF

		mergeit -p {output[0]} > {log} 2>&1
		"""


rule make_bed_Papac2021:
	input:
		multiext("data/Papac2021/tmp", ".bed", ".bim", ".fam"),
	output:
		multiext("data/Papac2021/Papac2021", ".bed", ".bim", ".fam"),
	params:
		prefix_in = lambda w, input: input[0].replace(".bed", ""),
		prefix_out = lambda w, output: output[0].replace(".bed", ""),
	conda: "../envs/eigensoft.yaml"
	shell:
		"""
		plink --bfile {params.prefix_in} \
		--make-bed --allow-no-sex --out {params.prefix_out}
		"""


# rule convert_plink2sgkit_Patterson2022:
# 	input:
# 		multiext("data/Patterson2022/Patterson2022", ".bed", ".bim", ".fam"),
# 		metadata_S3 = "data/Patterson2022/Patterson2022_TableS3.tsv",
# 		metadata_S5 = "data/Patterson2022/Patterson2022_TableS5_corrected.tsv",
# 		modern_anno = "data/v50.0_1240k_public_modern_selection.anno",
# 		bval_table = 'data/Murphy2021_Bvalues_compiled.bmap.txt',
# 		rval_table = 'data/Bherer2017_Refined_EUR_genetic_map_sexavg.rmap.txt',
# 	output:
# 		zarr_store = directory("data/Patterson2022/Patterson2022.zarr"),
# 		meta_out = "data/Patterson2022/Patterson2022_metadata.tsv",
# 		out_ind = "data/Patterson2022/Patterson2022.ind",
# 	params:
# 		path = lambda w, input: input[0].replace(".bed", ""),
# 	conda: "../envs/py-env.yaml"
# 	script:
# 		"../scripts/convert_plink2sgkit_Patterson2022.py"
