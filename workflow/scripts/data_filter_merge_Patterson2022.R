library(tidyverse)

samples_list <- read_lines(snakemake@input[['sample_list']])
modern <- read_lines(snakemake@input[['aadr_modern']])

aadr <- read_table(snakemake@input[['aadr_ind']],
				   col_names = c("sample", "sex", "pop"), col_types = 'ccc')
aadr %>% 
	mutate(pop = ifelse((aadr$sample %in% samples_list) | (aadr$sample %in% modern), pop, "Ignore")) %>%
	write_tsv(snakemake@output[[1]], col_names = F)

brit <- read_table(snakemake@input[['brit_ind']],
				   col_names = c("sample", "sex", "pop"), col_types = 'ccc')
brit %>%
	mutate(pop = ifelse(brit$sample %in% samples_list, pop, "Ignore")) %>%
	write_tsv(snakemake@output[[2]], col_names = F)

fowler <- read_table(snakemake@input[['hazleton_ind']],
				   col_names = c("sample", "sex", "pop"), col_types = 'ccc')
fowler %>%
	mutate(pop = ifelse(fowler$sample %in% samples_list, pop, "Ignore")) %>%
	write_tsv(snakemake@output[[3]], col_names = F)
