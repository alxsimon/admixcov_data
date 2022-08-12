import pandas as pd

bmap_files = snakemake.input['bmap_files']
rmap_files = snakemake.input['rmap_files']
chr_len_file = snakemake.input['chr_len']
bmap_out = snakemake.output['bmap_out']
rmap_out = snakemake.output['rmap_out']

df = pd.DataFrame()
for i, file in enumerate(bmap_files):
	new = pd.read_csv(file, sep='\s+', header=None, names=['bval', 'len'])
	new['start'] = new.len.cumsum() - new.len
	new['end'] = new.len.cumsum()
	new['chr'] = f'chr{i+1}'
	new = new[['chr', 'start', 'end', 'bval', 'len']]
	df = pd.concat([df, new], axis=0)
df.to_csv(bmap_out, sep='\t', index=False)


chr_len = pd.read_csv(chr_len_file, sep='\t', header=None, names=['chr', 'len'])
df = pd.DataFrame()
for file in rmap_files:
	new = pd.read_csv(file, sep='\t')
	new['start'] = new.pos - 1
	new['end'] = pd.concat([new.start[1:], chr_len.len[chr_len.chr == new.chr[0]]], ignore_index=True)
	df = pd.concat([df, new], axis=0)
	df = df[['chr', 'pos', 'start', 'end', 'rate', 'cM']]
df.to_csv(rmap_out, sep='\t', index=False)