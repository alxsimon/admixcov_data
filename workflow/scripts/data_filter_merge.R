library(tidyverse)

aadr <- read_table("data/v50.0_1240k_public.ind",
				   col_names = c("sample", "sex", "pop"), col_types = 'ccc')
TS3 <- read_tsv(snakemake@input[[1]], col_types = cols())

aadr_filt <- aadr %>% 
	mutate(pop = ifelse(aadr$sample %in% TS3$`Version ID`, pop, "Ignore"))

write_tsv(aadr_filt, "data/v50.0_1240k_public_filtered.ind", col_names = F)

writeLines(
	c(
		"geno1: data/brit.geno",
		"snp1: data/brit.snp",
		"ind1: data/brit.ind",
		"geno2: data/v50.0_1240k_public.geno",
		"snp2: data/v50.0_1240k_public.snp",
		"ind2: data/v50.0_1240k_public_filtered.ind",
		"outputformat: PACKEDPED",
		"genooutfilename: data/tmp_Patterson2022.bed",
		"snpoutfilename: data/tmp_Patterson2022.bim",
		"indoutfilename: data/tmp_Patterson2022.fam",
		"allowdups: NO"
	),
	"data/mergit_Patterson2022.par"
)
