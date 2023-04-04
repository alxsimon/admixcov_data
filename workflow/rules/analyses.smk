rule analysis_patterson:
	input:
		zarr = 'data/Patterson2022/Patterson2022.zarr',
	output:
		report = 'results/Patterson2022/analyses_info_Patterson2022.txt',
		fig = 'results/Patterson2022/fig_Patterson2022_main.pdf',
		fig_bins_G = 'results/Patterson2022/fig_Patterson2022_bins_G.pdf',
		fig_bins_var = 'results/Patterson2022/fig_Patterson2022_bins_var.pdf',
	conda:
		"../envs/py-env.yaml"
	script:
		'../scripts/analysis_patterson.py'

rule analysis_papac:
	input:
		zarr = 'data/Papac2021/Papac2021.zarr',
	output:
		report = 'results/Papac2021/analyses_info_Papac2021.txt',
		fig = 'results/Papac2021/fig_Papac2021_main.pdf',
		fig_bins_G = 'results/Papac2021/fig_Papac2021_bins_G.pdf',
		fig_bins_var = 'results/Papac2021/fig_Papac2021_bins_var.pdf',
	conda:
		"../envs/py-env.yaml"
	script:
		'../scripts/analysis_papac.py'