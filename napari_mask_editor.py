from pathlib import Path
import pandas as pd
import shutil
import tifffile
import matplotlib.pylab as plt
import matplotlib as mpl
mpl.rcParams["figure.dpi"] = 300
import numpy as np
import napari
from skimage import data
import xarray as xr
import re
import czifile

#%%
# LOAD DATASET
path_project = Path(r"Z:\0-Projects and Experiments\GG - ClickChemistry")


# list_str_path_generated_masks = list(map(str,list(path_project.rglob("*_mask_toxo.tiff"))))

list_str_path_generated_masks = list(map(str,list(path_project.rglob("*_cellpose.tiff"))))

list_czi_files = [str(p) for p in list(path_project.rglob("*.czi"))]


print(f"Number of czi files found: {len(list_czi_files)}")

index_start = 1# remember lists are zero index 
index_end =  3 # up to but not including

# ITERATE THROUGH ALL THE IMAGES
for path_mask in list_str_path_generated_masks[index_start -1 : index_end]:
    pass

    # CREATE AND MODIFY VIEWER
    viewer = napari.Viewer(show=False)
    
    # @viewer.bind_key("x", overwrite=True)
    # def exit(viewer):
    #     import napari
    #     from qtpy.QtCore import QTimer
    #     # with napari.gui_qt() as app:
    #         # viewer = napari.Viewer()
    #     time_in_msec = 1000
    #     QTimer().singleShot(time_in_msec, app.quit)
    #     # viewer.close()
    #     print("end")
            
    # @viewer.bind_key("s", overwrite=True)
    # def save_image(viewer):
    #     print("saving")
        
    # @viewer.bind_key("n", overwrite=True)
    # def next_image(viewer):
    #     print("next image")    
    
    # LOAD IMAGES
    path_mask = Path(path_mask)
    base_name = path_mask.stem.split("_")[0]
    path_czifile = list(filter(re.compile(f".*{base_name}.*").search, list_czi_files))[0]
    
    im_czi_data = czifile.imread(path_czifile).squeeze()
    im_czi_data = xr.DataArray(im_czi_data, dims=('c','x','y'),
                               coords={'c': ['toxo', 'dic', 'inosine', 'dapi']})
    

    #change to np.uint16 so we have more values for rois 2^16
    mask = tifffile.imread(path_mask).astype(np.uint16)
    
    # POPULATE VIEWER
    # layer_intensity = viewer.add_image(im_czi_data.sel(c='toxo'), name=Path(path_czifile).name)
    layer_intensity = viewer.add_image(im_czi_data.sel(c='inosine'), name=Path(path_czifile).name)

    layer_intensity.contrast_limits=(0,50)
    
    layer_mask = viewer.add_labels(mask, 
                                   name=path_mask.name,
                                   # color={1:'cyan'}
                                   )
    layer_mask.opacity = 0.4
    viewer.show(block=True)
    
    # # SAVE MASKS
    # path_im_output = path_mask / path_mask.name
    print(path_mask)
    tifffile.imwrite(path_mask, layer_mask.data)
    

    

    
    
    
    
    