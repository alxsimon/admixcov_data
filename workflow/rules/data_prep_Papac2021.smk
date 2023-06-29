
rule data_filter_Papac2021:
	input: 
		multiext("data/v50.0_1240k_public", ".geno", ".snp", ".ind"),
		"data/v50.0_1240k_public_modern_selection.anno",
		"data/Papac2021/Papac2021_TableS5_ancient_published.tsv",
		"data/Papac2021/Papac2021_TableS4_Bohemia_samples_info.tsv",
	output:
		"data/v50.0_1240k_public_filtered_Papac2021.ind",
	conda: "../envs/r-env.yaml"
	script:
		"../scripts/data_filter_merge_Papac2021.R"


rule filter_AADR_Papac2021:
	input:
		# multiext("data/Papac2021/Papac2021_forUpload", ".geno", ".modif.snp", ".ind"),
		multiext("data/v50.0_1240k_public", ".geno", ".snp", "_filtered_Papac2021.ind"),
	output:
		"data/Papac2021/convertf_bed_filter.par",
		temp(multiext("data/Papac2021/tmp", ".bed", ".bim", ".fam")),
	log:
		"logs/filter_AADR_Papac2021.log"
	conda: "../envs/eigensoft.yaml"
	shell:
		"""
		cat << EOF > {output[0]}
genotypename: {input[0]}
snpname: {input[1]}
indivname: {input[2]}
outputformat: PACKEDPED
genotypeoutname: {output[1]}
snpoutname: {output[2]}
indivoutname: {output[3]}
hashcheck: NO
EOF

		convertf -p {output[0]} > {log} 2>&1
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


rule convert_plink2sgkit_Papac2021:
	input:
		multiext("data/Papac2021/Papac2021", ".bed", ".bim", ".fam"),
		metadata_TS4 = "data/Papac2021/Papac2021_TableS4_Bohemia_samples_info.tsv",
		metadata_TS5 = "data/Papac2021/Papac2021_TableS5_ancient_published.tsv",
		metadata_TS9 = "data/Papac2021/Papac2021_TableS9.tsv",
		name_mapping = "data/Papac2021/Papac2021_sample_names_mapping.tsv",
		anno = "data/v50.0_1240k_public.anno",
		bval_table = 'data/Murphy2021_Bvalues_compiled.bmap.txt',
		rval_table = 'data/Bherer2017_Refined_EUR_genetic_map_sexavg.rmap.txt',
	output:
		zarr_store = directory("data/Papac2021/Papac2021.zarr"),
	params:
		path = lambda w, input: input[0].replace(".bed", ""),
	conda: "../envs/py-env.yaml"
	script:
		"../scripts/convert_plink2sgkit_Papac2021.py"
