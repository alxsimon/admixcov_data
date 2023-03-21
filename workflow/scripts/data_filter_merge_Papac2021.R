library(tidyverse)

# reformat SNP ids to match AADR
# snps_aadr = read_table(
# 	snakemake@input[[2]],
# 	col_names=c("SNP_ID", "Chromosome_Num", "Genetic_Position", "Physical_Position", "ref", "alt"),
# 	col_types="cidicc"
# ) |> unite('chr_pos', Chromosome_Num, Physical_Position, remove=FALSE)
# snps_papac = read_table(
# 	snakemake@input[[5]],
# 	col_names=c("chr_pos", "Chromosome_Num", "Genetic_Position", "Physical_Position", "ref", "alt"),
# 	col_types="cidicc"
# )
# snps_papac['SNP_ID'] = snps_aadr['SNP_ID'][match(snps_papac['chr_pos'], snps_aadr['chr_pos'])]
# snps_papac['Genetic_Position'] = snps_aadr['Genetic_Position'][match(snps_papac['chr_pos'], snps_aadr['chr_pos'])]
# snps_papac |> select(SNP_ID, Chromosome_Num, Genetic_Position, Physical_Position, ref, alt) |>
# 	write_tsv(snakemake@output[[1]], col_names=F)

# Papac2021 individuals are already in AADR v50.0

# Filter AADR individuals
aadr = read_table(snakemake@input[[3]],
	col_names=c("sample", "sex", "pop"), col_types='ccc')
# modern_papac = read_table(snakemake@input[[7]],
# 		col_names=c("sample", "sex", "pop"), col_types='ccc') |>
# 	filter(sample %in% aadr$sample)
modern = read_tsv(snakemake@input[[4]])
# Only take the references we are interested in
TS5 = read_tsv(snakemake@input[[5]]) |>
	filter(`Instance ID` %in% aadr$sample) |>
	filter(`Group_ID` %in% c(
		'Anatolia_Neolithic',
		'WHG (Loschbour)', 'WHG (Koros_HG)', 'WHG (Germany_Meso)',
		'Yamnaya_Samara', 'Yamnaya_Kalmykia'
	))
# Get Bohemian individuals
TS4 = read_tsv(snakemake@input[[6]]) |>
	filter(`Sample` %in% aadr$sample)

aadr |>
	mutate(pop = ifelse(
			(aadr$sample %in% modern$`Version ID`)
			| (aadr$sample %in% TS5$`Instance ID`)
			| (aadr$sample %in% TS4$Sample),
			pop, "Ignore"
		)) |>
	write_tsv(snakemake@output[[1]], col_names = F)
