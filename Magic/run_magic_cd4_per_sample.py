import magic
import os
import pandas as pd

rootdir = "data/th_cells_per_sample/samples"
outdir = "data/th_cells_per_sample/magic_output/"

for subdir, dirs, files in os.walk(rootdir):
    # loop over expression files
    for dir in dirs:
        scdata = magic.mg.SCData.from_10x(os.path.expanduser(os.path.join(subdir, dir)), 
                                          normalize=False)
        # run magic    
        scdata.run_magic(n_pca_components=20, random_pca=True, t=4, k=9, 
                 ka=3, epsilon=1, rescale_percent=99)
        
        # Take mean and add to eqtl_samples dataframe
        mean = scdata.magic.data.mean()
        mean.name = dir
        
        # write imputed expression file per sample
        scdata.magic.data.transpose().to_csv(os.path.join(outdir, dir + ".tsv"), header=True, index=True, sep='\t', mode='w')

