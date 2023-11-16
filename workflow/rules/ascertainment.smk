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
		# """
		# rg -N -z \
		# -f <(awk -v chr="{wildcards.chr}" '$1==chr {{print "^"$1"\\\s"$2"\\\s"}}' {input.used_snps}) \
		# {input.rohland_chr} \
		# > {output.out}
		# """


