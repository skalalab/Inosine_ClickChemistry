from pathlib import Path

import matplotlib.pylab as plt
import matplotlib as mpl
mpl.rcParams['figure.dpi'] = 300

import tifffile
import czifile

import numpy as np
import re

import cell_analysis_tools as cat
from cell_analysis_tools.visualization import compare_images
from cell_analysis_tools.image_processing import kmeans_threshold
#%%

path_project = Path(r"Z:\0-Projects and Experiments\GG - ClickChemistry")
list_czi_files = list(path_project.rglob("*.czi"))


for path_czi in list_czi_files[:]:
    pass
    
    base_name = path_czi.stem
    list_folders_in_dir = [str(p) for p in path_czi.parent.glob("*") if p.is_dir()]

    im = czifile.imread(path_czi).squeeze()
    
    im_toxo_mask = im[0,...]
    toxo_mask = kmeans_threshold(im_toxo_mask, 
                                 k=4,
                                 n_brightest_clusters=2)
    
    compare_images('original', im_toxo_mask, 
                   'mask',toxo_mask,
                   suptitle=f"{path_czi.name}")
    
    
    path_output = list(filter(re.compile(f".*{base_name}.*").search, list_folders_in_dir))[0]
    # save mask
    
    tifffile.imwrite(Path(path_output) / f"{base_name}_mask_toxo.tiff", toxo_mask)
