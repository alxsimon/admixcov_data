
rule smartpca:
	input:
		multiext("data/Patterson2022", ".bed", ".bim", ".ind"),
	output:
		multiext("results/smartpca/full_pca", ".evec", ".eval", ".snpweight"),
		"results/smartpca/full_pca.params",
	log:
		"logs/smartpca_full_pca.log"
	conda: "../envs/analyses.yaml"
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

rule smartpca_ref:
	input:
		multiext("data/Patterson2022", ".bed", ".bim", ".ind"),
	output:
		multiext("results/smartpca/ref_pca", ".evec", ".eval", ".snpweight"),
		"results/smartpca/ref_pca.params",
		"results/smartpca/ref_pca_used_groups.txt",
	log:
		"logs/smartpca_ref_pca.log"
	conda: "../envs/analyses.yaml"
	threads:
		workflow.cores
	shell:
		"""
		printf 'WHGA\\nBalkan_N\\nOldSteppe\\n' > {output[4]}

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
poplistname: {output[4]}
fsthiprecision: YES
inbreed: YES
numthreads: {threads}
EOF

		smartpca -p {output[3]} > {log} 2>&1
		"""