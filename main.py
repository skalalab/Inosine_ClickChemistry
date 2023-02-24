from pathlib import Path

import matplotlib.pylab as plt
import matplotlib as mpl
mpl.rcParams['figure.dpi'] = 300

import tifffile
import czifile

import numpy as np

import cell_analysis_tools as cat
from cell_analysis_tools.visualization import compare_images
from cell_analysis_tools.image_processing import kmeans_threshold
#%%

path_project = Path(r"Z:\0-Projects and Experiments\GG - ClickChemistry")

list_czi_files = list(path_project.rglob("*.czi"))


for path_czi in list_czi_files[:10]:
    pass

    im = czifile.imread(path_czi).squeeze()
    
    im_toxo_mask = im[0,...]
    percentile_threshold = np.percentile(im_toxo_mask,99)
    toxo_mask = im_toxo_mask > (im_toxo_mask.max()-3)

    toxo_mask = kmeans_threshold(im_toxo_mask, 
                                 k=4, 
                                 n_brightest_clusters=1)
    
    compare_images('original', im_toxo_mask, 
                   'mask',toxo_mask,
                   suptitle=f"{path_czi.name}")
    
    
    # for data in im:
    #     plt.title(f"{path_czi.name}")
    #     plt.imshow(data)
    #     plt.show()
