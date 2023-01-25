
rule smartpca_Patterson2022:
	input:
		multiext("data/Patterson2022/Patterson2022", ".bed", ".bim", ".ind"),
	output:
		multiext("results/smartpca_Patterson2022/full_pca", ".evec", ".eval", ".snpweight"),
		"results/smartpca_Patterson2022/full_pca.params",
	log:
		"logs/smartpca_Patterson2022_full_pca.log"
	conda: "../envs/eigensoft.yaml"
	threads:
		workflow.cores
	shell:
		"""
		cat << EOF > {output[3]}
genotypename: {input[0]}
snpname: {input[1]}
indivname: {input[2]}
evecoutname: {output[0]}
evaloutname: {output[1]}
snpweightoutname: {output[2]}
numoutevec: 10
numoutlieriter: 0
fsthiprecision: YES
inbreed: YES
numthreads: {threads}
EOF

		smartpca -p {output[3]} > {log} 2>&1
		"""

rule smartpca_Patterson2022_proj:
	input:
		multiext("data/Patterson2022/Patterson2022", ".bed", ".bim", ".ind"),
		"data/Patterson2022/pca_projection_modern_groups.txt",
	output:
		multiext("results/smartpca_Patterson2022/pca_proj", ".evec", ".eval", ".snpweight"),
		"results/smartpca_Patterson2022/pca_proj.params",
		# "results/smartpca/pca_proj_used_groups.txt",
	log:
		"logs/smartpca_Patterson2022_proj.log"
	conda: "../envs/eigensoft.yaml"
	threads:
		workflow.cores
	shell:
		"""
		cat << EOF > {output[3]}
genotypename: {input[0]}
snpname: {input[1]}
indivname: {input[2]}
evecoutname: {output[0]}
evaloutname: {output[1]}
snpweightoutname: {output[2]}
numoutevec: 10
numoutlieriter: 0
lsqproject: YES
poplistname: {input[3]}
fsthiprecision: YES
inbreed: YES
numthreads: {threads}
EOF

		smartpca -p {output[3]} > {log} 2>&1
		"""