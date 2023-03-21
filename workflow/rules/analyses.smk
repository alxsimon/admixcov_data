rule analysis_patterson:
	input:
		zarr = 'data/Patterson2022/Patterson2022.zarr',
	output:
		fig = 'results/Patterson2022/main_figure_Patterson2022.pdf',
	conda:
		"../envs/py-env.yaml"
	script:
		'../scripts/analysis_patterson.py'

rule analysis_papac:
	input:
		zarr = 'data/Papac2021/Papac2021.zarr',
	output:
		fig = 'results/Papac2021/main_figure_Papac2021.pdf',
	conda:
		"../envs/py-env.yaml"
	script:
		'../scripts/analysis_papac.py'