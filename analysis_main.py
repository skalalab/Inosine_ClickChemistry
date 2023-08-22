from pathlib import Path

import matplotlib.pylab as plt
import matplotlib as mpl
mpl.rcParams['figure.dpi'] = 300

import tifffile
import czifile

from natsort import natsorted

import numpy as np
import re

from skimage.morphology import label

from skimage.filters import threshold_multiotsu

import cell_analysis_tools as cat
from cell_analysis_tools.visualization import compare_images
from cell_analysis_tools.image_processing import kmeans_threshold

from tqdm import tqdm
import pandas as pd
#%%

path_project = Path(r"Z:\0-Projects and Experiments\GG - ClickChemistry")
list_czi_files = list(path_project.rglob("*.czi"))


df_toxo_inosine = pd.DataFrame()
df_cell_inosine = pd.DataFrame()


for idx, path_czi in tqdm(enumerate(list_czi_files[:45])): # threshold_multiosu slow on 11:12 
    pass
    base_name = path_czi.stem
    list_folders_in_dir = [str(p) for p in path_czi.parent.glob("*") if p.is_dir()]

    im = czifile.imread(path_czi).squeeze()
    
    bool_show_images = False
    
    # print show different channels and their indices
    if bool_show_images:
        fig, ax = plt.subplots(1,4, figsize=(10,3))
        fig.suptitle(f"{path_czi}")
        for idx in range(len(im)):
            pass
            ax[idx].set_title(f"index: {idx}")
            ax[idx].set_axis_off()
            ax[idx].imshow(im[idx,...])
        plt.show()
            
            # plt.title(f"index: {idx}")
            # plt.imshow(im[idx,...])
            # plt.show()
    
    # Channels:
    ch_toxo = 0 # red 
    ch_DIC = 1 # 
    ch_inosine = 2 # green
    ch_DAPI = 3 # DIC HFF
    
    im_toxo = im[ch_toxo,...]
    
    ### NOTE: These lines of code were how the initial toxo masks were created 
    ### and then were revised by Gina
    # toxo_mask = kmeans_threshold(im_toxo, 
    #                              k=4,
    #                              n_brightest_clusters=2)
    
    # This creates initial toxo mask that is revised by Gina  
    ##### tifffile.imwrite(Path(path_output) / f"{base_name}_mask_toxo.tiff", toxo_mask)
    
    
    #### inosine 
    im_inosine = im[ch_inosine,...]
    
    # thresh_inosine = threshold_multiotsu(im_inosine, 4)
    
    # mask_inosine = im_inosine > thresh_inosine[-1:]
    # if bool_show_images:
    #     compare_images('original inosine', im_inosine,
    #                     'mask', mask_inosine,
    #                     suptitle=f"\n{path_czi}")
        
    
    # using multiotsu
    # thresh_inosine = threshold_multiotsu(im_inosine, 7)
    # mask_inosine = im_inosine > thresh_inosine[-1:]
    
    
    # set specific k and keep prameters for images
    dict_parameters = {
        # 'Snap-483' : {'k': 7, 'keep': 2},
        # 'Snap-482' : {'k': 6, 'keep': 1},
        # 'Snap-494' : {'k': 7, 'keep': 1},
        # 'Snap-497' : {'k': 7, 'keep': 2},
        'Snap-679' : {'k': 7, 'keep': 1},
        }
    
    if base_name in dict_parameters:
        k = dict_parameters[base_name]['k']
        keep = dict_parameters[base_name]['keep']
    else: # general paraks for rest of image 
    # using kmeans 
        k = 7
        keep = 2
        

    # CREATE INOSINE MASK 
    mask_inosine = kmeans_threshold(im_inosine, k=k, n_brightest_clusters=keep)
    
    if bool_show_images:
        compare_images('original inosine', im_inosine,
                        'mask', mask_inosine,
                        suptitle=f"{idx} | \n{path_czi}")

    # LOAD TOXO MASK
    path_output = list(filter(re.compile(f".*{base_name}.*").search, list_folders_in_dir))[0]
    # save mask
    
    mask_toxo = tifffile.imread(Path(path_output) / f"{base_name}_mask_toxo.tiff")
    mask_toxo = label(mask_toxo)
    
    if bool_show_images:
        compare_images('im_toxo', im_toxo,
                        'toxo mask', mask_toxo, 
                        suptitle=f"{path_czi.name}")

    # LOAD WHOLE CELL MASK
    mask_cell = tifffile.imread(Path(path_output).parent / f"{base_name}_cellpose.tiff")
    mask_cell = label(mask_cell)
    
    
    _, host_cell, toxo_strain = path_czi.parent.parent.name.split(" ")

    
    ############### COMPUTE WHOLE CELL INOSINE QUANTIFICATION

    # create cell mask without toxo
    mask_cell_no_toxo = mask_cell * np.invert(mask_toxo > 0)
    
    # get unique labels for remaining regions
    list_cell_labels = list(np.unique(mask_cell_no_toxo))
    list_cell_labels.remove(0)
    
    for label_cell in list_cell_labels:
        pass
    
        mask_single_cell = mask_cell == label_cell
        
        area_cell = np.sum(mask_single_cell)
        mask_intra_cell_inosine = mask_single_cell * mask_inosine
        area_intra_cell_inosine = np.sum(mask_intra_cell_inosine)
        
        percent_intra_cell_inosine = area_intra_cell_inosine / area_cell
        intensity_sum_intra_cell_inosine = np.sum(area_intra_cell_inosine * im_inosine)
        intensity_mean_intra_cell_inosine =  intensity_sum_intra_cell_inosine/ area_cell
        
        
        dict_properties = {
            'file_name' : path_czi.name,
            'label_cell' : label_cell,
            'path_file' : str(path_czi),
            'replicate' : int(path_czi.parent.name[:1]),
            'Host Cell' : host_cell,
            'Toxoplasma Strain' : toxo_strain,
            'area_cell' : area_cell,
            'area_intra_cell_inosine' : area_intra_cell_inosine,
            'intensity_sum_intra_cell_inosine': intensity_sum_intra_cell_inosine,
            'intensity_mean_intra_toxo_inosine' : intensity_mean_intra_cell_inosine,
            'percent_intra_cell_inosine' : area_intra_cell_inosine / area_cell,
            }
        
        df_single_cell_roi = pd.DataFrame(dict_properties, index=[0])
        
        df_cell_inosine = pd.concat([df_cell_inosine, df_single_cell_roi], ignore_index=True)


    ############### COMPUTE INTRA TOXO INOSINE
    list_toxo_labels = list(np.unique(mask_toxo))
    list_toxo_labels.remove(0)
    
    for label_toxo in list_toxo_labels:
        pass
    
        # populate dataframe
        
        # generate overlap mask
        mask_single_toxo = mask_toxo == label_toxo
        mask_intra_toxo_inosine = mask_single_toxo * mask_inosine
        
        # compute area overlap
        area_toxo = np.sum(mask_single_toxo)
        area_intra_toxo_inosine = np.sum(mask_intra_toxo_inosine)
        percent_intra_toxo_inosine = area_intra_toxo_inosine / area_toxo
        
        # integrated intensity of insoine / area of toxo
        intensity_intra_toxo_inosine = mask_intra_toxo_inosine * im_inosine
        intensity_mean_intra_toxo_inosine = np.sum(intensity_intra_toxo_inosine) / area_toxo

        ## visualize overlap
        if bool_show_images:
            compare_images('single toxo mask', mask_single_toxo, 
                            'mask_intra_inosine', mask_intra_toxo_inosine,
                            suptitle=f"{path_czi.name} \ntoxo_label={label_toxo} | percent_intra_inosine: {percent_intra_toxo_inosine:.3f}")
        
        dict_properties = {
            'file_name' : path_czi.name,
            'label_toxo' : label_toxo,
            'path_file' : str(path_czi),
            'replicate' : int(path_czi.parent.name[:1]),
            'Host Cell' : host_cell,
            'Toxoplasma Strain' : toxo_strain,
            'area_inosine' : area_intra_toxo_inosine,
            'area_toxo' : area_toxo,
            'intensity_intra_toxo_inosine': np.sum(intensity_intra_toxo_inosine),
            'percent_intra_toxo_inosine' : percent_intra_toxo_inosine,
            'intensity_mean_intra_toxo_inosine' : intensity_mean_intra_toxo_inosine,
            }
        
        df_single_roi = pd.DataFrame(dict_properties, index=[0])
        
        df_toxo_inosine = pd.concat([df_toxo_inosine, df_single_roi], ignore_index=True)
        
    # save intra toxo inosine df
    filename = f"click_chemistry_intra_toxo_inosine_features.csv"
    df_toxo_inosine.to_csv(path_project / filename)
    
    # save intra cell inosine df
    filename = f"click_chemistry_intra_cell_inosine_features.csv"
    df_cell_inosine.to_csv(path_project / filename)
        
        

#%% plots


import seaborn as sns

df_toxo_inosine['conditions'] = df_toxo_inosine['Host Cell'] + "_"+ df_toxo_inosine['Toxoplasma Strain']

ax = sns.swarmplot(df_toxo_inosine, x='conditions', y='percent_inosine_in_toxo', hue='replicate')
plt.suptitle(f"inosine k: {k} | keep: {keep}")
sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))










