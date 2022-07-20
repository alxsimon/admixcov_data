library(tidyverse)

samples_list <- read_lines(snakemake@input[[1]])

aadr <- read_table(snakemake@input[[2]],
				   col_names = c("sample", "sex", "pop"), col_types = 'ccc')
aadr %>% 
	mutate(pop = ifelse(aadr$sample %in% samples_list, pop, "Ignore")) %>%
	write_tsv("data/v50.0_1240k_public_filtered.ind", col_names = F)

brit <- read_table(snakemake@input[[3]],
				   col_names = c("sample", "sex", "pop"), col_types = 'ccc')
brit %>%
	mutate(pop = ifelse(brit$sample %in% samples_list, pop, "Ignore")) %>%
	write_tsv("data/brit_filtered.ind", col_names = F)

fowler <- read_table(snakemake@input[[4]],
				   col_names = c("sample", "sex", "pop"), col_types = 'ccc')
fowler %>%
	mutate(pop = ifelse(fowler$sample %in% samples_list, pop, "Ignore")) %>%
	write_tsv("data/Hazleton_v51_filtered.ind", col_names = F)

writeLines(
	c(
		"geno1: data/brit.geno",
		"snp1: data/brit.snp",
		"ind1: data/brit_filtered.ind",
		"geno2: data/Hazleton_v51.geno",
		"snp2: data/Hazleton_v51.snp",
		"ind2: data/Hazleton_v51_filtered.ind",
		"outputformat: EIGENSTRAT",
		"genooutfilename: data/tmp_1.geno",
		"snpoutfilename: data/tmp_1.snp",
		"indoutfilename: data/tmp_1.ind",
		"allowdups: NO",
		"hashcheck: NO"
	),
	"data/mergit_intermediate.par"
)

writeLines(
	c(
		"geno1: data/tmp_1.geno",
		"snp1: data/tmp_1.snp",
		"ind1: data/tmp_1.ind",
		"geno2: data/v50.0_1240k_public.geno",
		"snp2: data/v50.0_1240k_public.snp",
		"ind2: data/v50.0_1240k_public_filtered.ind",
		"outputformat: PACKEDPED",
		"genooutfilename: data/tmp_2.bed",
		"snpoutfilename: data/tmp_2.bim",
		"indoutfilename: data/tmp_2.fam",
		"allowdups: NO"
	),
	"data/mergit_Patterson2022.par"
)
