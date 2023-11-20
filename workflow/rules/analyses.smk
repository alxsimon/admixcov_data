rule analysis_patterson:
	input:
		zarr = 'data/Patterson2022/Patterson2022.zarr',
	output:
		report = 'results/Patterson2022/analyses_info_Patterson2022.txt',
		fig_data = 'results/Patterson2022/fig_data_Patterson2022_main.pickle',
		fig_complete_G = 'results/Patterson2022/fig_Patterson2022_complete_G.pdf',
		fig_bins_G = 'results/Patterson2022/fig_Patterson2022_bins_G.pdf',
		fig_bins_var = 'results/Patterson2022/fig_Patterson2022_bins_var.pdf',
		fig_bins_totvar = 'results/Patterson2022/fig_Patterson2022_bins_totvar.pdf',
		matrix_data = 'results/Patterson2022/matrix_Patterson2022.pickle',
		fig_data_bval_bins = 'results/Patterson2022/fig_data_Patterson_bval_bins.pickle',
		used_snps = 'results/Patterson2022/used_snps_Patterson2022.tsv',
	conda:
		"../envs/py-env.yaml"
	resources:
		mem = "10G"
	script:
		'../scripts/analysis_patterson.py'

rule analysis_papac:
	input:
		zarr = 'data/Papac2021/Papac2021.zarr',
	output:
		report = 'results/Papac2021/analyses_info_Papac2021.txt',
		fig_data = 'results/Papac2021/fig_data_Papac2021_main.pickle',
		fig_complete_G = 'results/Papac2021/fig_Papac2021_complete_G.pdf',
		fig_bins_G = 'results/Papac2021/fig_Papac2021_bins_G.pdf',
		fig_bins_var = 'results/Papac2021/fig_Papac2021_bins_var.pdf',
		fig_bins_totvar = 'results/Papac2021/fig_Papac2021_bins_totvar.pdf',
		matrix_data = 'results/Papac2021/matrix_Papac2021.pickle',
		used_snps = 'results/Papac2021/used_snps_Papac2021.tsv',
	conda:
		"../envs/py-env.yaml"
	resources:
		mem = "10G"
	script:
		'../scripts/analysis_papac.py'

rule main_figures:
	input:
		matrix_uk = 'results/Patterson2022/matrix_Patterson2022.pickle',
		matrix_bo = 'results/Papac2021/matrix_Papac2021.pickle',
		file_data_uk = 'results/Patterson2022/fig_data_Patterson2022_main.pickle',
		file_data_bo = 'results/Papac2021/fig_data_Papac2021_main.pickle',
		file_data_uk_bval_bins = 'results/Patterson2022/fig_data_Patterson_bval_bins.pickle',
	output:
		fig_data_matrix_ancestry = 'results/fig_data_matrix_ancestry.pdf',
		fig_data_covlines = 'results/fig_data_covlines.pdf',
		fig_data_uk_bval = 'results/fig_data_uk_bval.pdf',
	conda:
		"../envs/py-env.yaml"
	script:
		'../scripts/main_figures.py'

# split time intervals
rule analysis_patterson_split:
	input:
		zarr = 'data/Patterson2022/Patterson2022.zarr',
	output:
		report = 'results/Patterson2022_split/analyses_info_Patterson2022_split.txt',
		fig = 'results/Patterson2022_split/fig_Patterson2022_split_main.pdf',
		fig_bins_G = 'results/Patterson2022_split/fig_Patterson2022_split_bins_G.pdf',
		fig_bins_var = 'results/Patterson2022_split/fig_Patterson2022_split_bins_var.pdf',
		fig_bins_totvar = 'results/Patterson2022_split/fig_Patterson2022_split_bins_totvar.pdf',
		matrix_data = 'results/Patterson2022_split/matrix_Patterson2022_split.pickle',
	conda:
		"../envs/py-env.yaml"
	resources:
		mem = "10G"
	script:
		'../scripts/analysis_patterson_split.py'
