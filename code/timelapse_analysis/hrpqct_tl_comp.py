"""
Script Name: hrpqct_tl_comp.py
Description: This script:
                - Runs the timelapse calculations and stores the results and the images
"""
from utils.logger_setup import logger
import os
import csv
import numpy as np
import matplotlib.pyplot as plt
import scipy
import imageio
from skimage.morphology import flood_fill
from skimage import util 
from skimage.morphology import label
from skimage.morphology import remove_small_objects
from skimage.morphology import binary_dilation
from skimage.morphology import binary_erosion
from skimage.morphology import (disk)
from datetime import date
from datetime import datetime
import time
from utils import aim_reader_v1 as aim


def map2common_frame(
    gray1,gray2,seg1,seg2,cort1,cort2,trab1,trab2,
    header_gray1,header_gray2,header_seg1,header_seg2,
    header_cort1,header_cort2,header_trab1,header_trab2
):
    """
    Function for aligning all image components into a common grid, to ensure proper alignment
    """
    nps = [gray1,gray2,seg1,seg2,cort1,cort2,trab1,trab2]   # Numpy arrays of images and masks
    headers = [header_gray1,header_gray2,header_seg1,header_seg2,header_cort1,\
               header_cort2,header_trab1,header_trab2]   # Headers of images
    
    logger.info("       Calculations start now!")
    
    if header_gray1["Patient_ID"] != header_gray2["Patient_ID"]:
        logger.info("Error! The patient ID does not match between the two scans!")
    Patient_ID = header_gray1["Patient_ID"]
    # ===========================================================================
    # Finding the number of days between scans. This will swap scans if they
    # are not ordered correctly according to their dates
    n_days = date(int(header_gray2["Date"][2]),int(header_gray2["Date"][1]),int(header_gray2["Date"][0]))\
        - date(int(header_gray1["Date"][2]),int(header_gray1["Date"][1]),int(header_gray1["Date"][0]))
    n_days = n_days.days
    if n_days == 0: n_days = 1
    elif n_days < 0:
        logger.info("Scans are not ordered chronologically. They will be swapped!")
        gray1,header_gray1 = aim.aim2np(gray2)
        gray2,header_gray2 = aim.aim2np(gray1)
        seg1,header_seg1 = aim.aim2np(seg2)
        seg2,header_seg2 = aim.aim2np(seg1)
        cort1,header_cort1 = aim.aim2np(cort2)
        cort2,header_cort2 = aim.aim2np(cort1)
        trab1,header_trab1 = aim.aim2np(trab2)
        trab2,header_trab2 = aim.aim2np(trab1)
        
        nps = [gray1,gray2,seg1,seg2,cort1,cort2,trab1,trab2]
        headers = [header_gray1,header_gray2,header_seg1,header_seg2,header_cort1,\
                   header_cort2,header_trab1,header_trab2]
    # ===========================================================================
    # Determine the voxel size to be used for calculations
    voxel_size = 82 if int(header_gray1["XCT_gen"]) == 1 else  60.7
    # =============================================================================
    # Finding the smallest 3D arrays that fits all of the images and aligning them in x-y plane
    order = ["gray1","gray2","seg1","seg2","cort1","cort2","trab1","trab2"]
    order_row = order[:]
    order_col = order[:]
    order_z = order[:]
    # Getting 3D coordinte of each image in the global coordinate system
    rows = [int(headers[i]["coordinates"][1]) for i in range(len(headers))]
    columns = [int(headers[i]["coordinates"][0]) for i in range(len(headers))]
    zs = [int(headers[i]["coordinates"][2]) for i in range(len(headers))]
    
    # Finding the highest value of x, y, and z coordinates
    rows_max = [int(headers[i]["coordinates"][1]+headers[i]["dimensions"][1])\
                for i in range(len(headers))]
    columns_max = [int(headers[i]["coordinates"][0]+headers[i]["dimensions"][0])\
                   for i in range(len(headers))]
    zs_max = [int(headers[i]["coordinates"][2]+headers[i]["dimensions"][2])\
              for i in range(len(headers))]
        
    # Sort the images based on their locations in x, y, and z axes
    tuples = zip(*sorted(zip(rows, order_row)))
    rows, order_row = [ list(tuple) for tuple in  tuples]
    tuples = zip(*sorted(zip(columns, order_col)))
    columns, order_col = [ list(tuple) for tuple in  tuples]
    tuples = zip(*sorted(zip(zs, order_z)))
    zs, order_z = [ list(tuple) for tuple in  tuples]
    
    row_c = max(rows_max) - rows[0]
    col_c = max(columns_max) - columns[0]
    z_c = max(zs_max) - zs[0]
    
    rows[:] = [elem - rows[0] for elem in rows]
    columns[:] = [elem - columns[0] for elem in columns]
    zs[:] = [elem - zs[0] for elem in zs]
    
    # Create zero numpy arrays with the common size for all images
    gray1_common = np.zeros((row_c,col_c,z_c), dtype="float32")     # Float becuase of the Gaussian filter
    gray2_common = np.zeros((row_c,col_c,z_c), dtype="float32")     # Float becuase of the Gaussian filter
    seg1_common = np.zeros((row_c,col_c,z_c), dtype="uint8")
    seg2_common = np.zeros((row_c,col_c,z_c), dtype="uint8")
    cort1_common = np.zeros((row_c,col_c,z_c), dtype="uint8")
    cort2_common = np.zeros((row_c,col_c,z_c), dtype="uint8")
    trab1_common = np.zeros((row_c,col_c,z_c), dtype="uint8")
    trab2_common = np.zeros((row_c,col_c,z_c), dtype="uint8")

    # Map the images to the common array
    gray1_common[rows[order_row.index("gray1")]:(rows[order_row.index("gray1")]+headers[0]["dimensions"][1])\
          ,columns[order_col.index("gray1")]:(columns[order_col.index("gray1")]+headers[0]["dimensions"][0])\
              ,zs[order_z.index("gray1")]:zs[order_z.index("gray1")]+headers[0]["dimensions"][2]] = np.int16(gray1)
        
    gray2_common[rows[order_row.index("gray2")]:(rows[order_row.index("gray2")]+headers[1]["dimensions"][1])\
          ,columns[order_col.index("gray2")]:(columns[order_col.index("gray2")]+headers[1]["dimensions"][0])\
              ,zs[order_z.index("gray2")]:zs[order_z.index("gray2")]+headers[1]["dimensions"][2]] = np.int16(gray2)
        
    seg1_common[rows[order_row.index("seg1")]:(rows[order_row.index("seg1")]+headers[2]["dimensions"][1])\
          ,columns[order_col.index("seg1")]:(columns[order_col.index("seg1")]+headers[2]["dimensions"][0])\
              ,zs[order_z.index("seg1")]:zs[order_z.index("seg1")]+headers[2]["dimensions"][2]] = np.uint(seg1)
        
    seg2_common[rows[order_row.index("seg2")]:(rows[order_row.index("seg2")]+headers[3]["dimensions"][1])\
          ,columns[order_col.index("seg2")]:(columns[order_col.index("seg2")]+headers[3]["dimensions"][0])\
              ,zs[order_z.index("seg2")]:zs[order_z.index("seg2")]+headers[3]["dimensions"][2]] = np.uint8(seg2)
        
    cort1_common[rows[order_row.index("cort1")]:(rows[order_row.index("cort1")]+headers[4]["dimensions"][1])\
          ,columns[order_col.index("cort1")]:(columns[order_col.index("cort1")]+headers[4]["dimensions"][0])\
              ,zs[order_z.index("cort1")]:zs[order_z.index("cort1")]+headers[4]["dimensions"][2]] = np.uint8(cort1)
        
    cort2_common[rows[order_row.index("cort2")]:(rows[order_row.index("cort2")]+headers[5]["dimensions"][1])\
          ,columns[order_col.index("cort2")]:(columns[order_col.index("cort2")]+headers[5]["dimensions"][0])\
              ,zs[order_z.index("cort2")]:zs[order_z.index("cort2")]+headers[5]["dimensions"][2]] = np.uint8(cort2)
        
    trab1_common[rows[order_row.index("trab1")]:(rows[order_row.index("trab1")]+headers[6]["dimensions"][1])\
          ,columns[order_col.index("trab1")]:(columns[order_col.index("trab1")]+headers[6]["dimensions"][0])\
              ,zs[order_z.index("trab1")]:zs[order_z.index("trab1")]+headers[6]["dimensions"][2]] = np.uint8(trab1)
        
    trab2_common[rows[order_row.index("trab2")]:(rows[order_row.index("trab2")]+headers[7]["dimensions"][1])\
          ,columns[order_col.index("trab2")]:(columns[order_col.index("trab2")]+headers[7]["dimensions"][0])\
              ,zs[order_z.index("trab2")]:zs[order_z.index("trab2")]+headers[7]["dimensions"][2]] = np.uint8(trab2)
    
    logger.info("All the image components are mapped to a common coordinate!")
    # =====================================================================
    # Masking the grayscale images based on the compartment to be analyzed
    cort1_common[np.where(trab1_common != 0)] = 0
    cort2_common[np.where(trab2_common != 0)] = 0
    
    # Making sure all masks are binary (0 and 1) instead of 126 or 127
    cort1_common[np.where(cort1_common != 0)] = 1
    cort2_common[np.where(cort2_common != 0)] = 1
    trab1_common[np.where(trab1_common != 0)] = 1 
    trab2_common[np.where(trab2_common != 0)] = 1    
    seg1_common[np.where(seg1_common != 0)] = 1    
    seg2_common[np.where(seg2_common != 0)] = 1
    
    return gray1_common,gray2_common,seg1_common,seg2_common,cort1_common,cort2_common,trab1_common,trab2_common,header_gray1,header_gray2

def gaussian_filter_3d(image,sigma=0.8):
    if sigma != 0:    
        filtered = scipy.ndimage.filters.gaussian_filter(image, sigma, order=0, output=None, mode='nearest', cval=0.0)
    elif sigma == 0:
        filtered = image
    return filtered

def gaussian_filter_3d_trunc(image,sigma=0.8):
    if sigma != 0:    
        filtered = scipy.ndimage.filters.gaussian_filter(image, sigma, order=0, output=None, mode='nearest', cval=0.0, truncate=1)
    elif sigma == 0:
        filtered = image
    return filtered

def gray_clean_up(gray1_common,gray2_common):
    """
    Function for removing voxels where one of the grayscale images is zero, otherwise the
    dges of partial cross-sections are usually counted as resorption due to aliasing
    """
    gray1_common[np.where(gray2_common == 0)] = 0
    gray2_common[np.where(gray1_common == 0)] = 0
    logger.info("Aliasing at the edges is fixed!")
    return gray1_common,gray2_common

def gray_subtraction(gray1_common,gray2_common,header_gray1,header_gray2,thresh = 225):
    """
    Function for subtracting the grayscale image to identify canidate formation/resorption voxels 
    """
    # ========================================================================
    # Converting gray native values to bone density mgHA/cm3
    gray1_common = gray1_common/(float(header_gray1["Mu_scaling"]))*header_gray1["slope"] + header_gray1["intercept"]
    gray2_common = gray2_common/(float(header_gray2["Mu_scaling"]))*header_gray2["slope"] + header_gray2["intercept"]
    histomorph = np.zeros_like(gray1_common)  
    histomorph = histomorph.astype("int16")
    tmp = np.subtract(gray1_common,gray2_common)
    
    histomorph[np.where(tmp > thresh)] = 500    # Resorption
    histomorph[np.where(tmp < -thresh)] = 300   # Formation
    return histomorph

def binary_subtraction(seg1_common,seg2_common):
    """
    Function for subracting binary images to identify formation/esorption voxels
    """
    seg1_common[np.where(seg1_common != 0)] = 1    
    seg2_common[np.where(seg2_common != 0)] = 1
    tmp = np.subtract(seg1_common,seg2_common)
    histomorph = np.zeros_like(seg1_common)  
    histomorph=histomorph.astype("int16")
    histomorph[np.where(tmp == 1)] = 500    # Resorption
    histomorph[np.where(tmp == 255)] = 300  # Formation
    return histomorph

def mask_images(
    gray1_common,gray2_common,seg1_common,seg2_common,
    cort1_common,cort2_common,trab1_common,trab2_common,\
         input_image,component,mask,cort_mask,cort_surface
): # total,cort,trab
    """
    Function for masking the desired bone compartment
    """
    # Mask the desired bone region
    cort1_common[np.where(trab1_common != 0)] = 0
    cort2_common[np.where(trab2_common != 0)] = 0
    cort1_common[np.where(cort1_common != 0)] = 1
    cort2_common[np.where(cort2_common != 0)] = 1
    trab1_common[np.where(trab1_common != 0)] = 1 
    trab2_common[np.where(trab2_common != 0)] = 1    
    seg1_common[np.where(seg1_common != 0)] = 1    
    seg2_common[np.where(seg2_common != 0)] = 1
    
    if component == "Total":
        mask1_common = cort1_common + trab1_common
        mask2_common = cort2_common + trab2_common
        mask1_common = util.invert(flood_fill(mask1_common, (1, 1, 1), 1)) + 2 + mask1_common
        mask2_common = util.invert(flood_fill(mask2_common, (1, 1, 1), 1)) + 2 + mask2_common
        mask1_common[np.where(mask1_common != 0)] = 1    
        mask2_common[np.where(mask2_common != 0)] = 1
        if input_image == "Grayscale":
            if mask == "same":
                struc_2d = {"disk": disk(6)}
                for i in range(np.shape(mask1_common)[2]):
                    mask1_common[:,:,i] = binary_dilation(mask1_common[:,:,i],struc_2d["disk"])
                gray1_common = (mask1_common) * gray1_common
                gray2_common = (mask1_common) * gray2_common
            elif mask == "different":
                gray1_common = (mask1_common) * gray1_common
                gray2_common = (mask2_common) * gray2_common
        elif input_image == "Binary":
            if mask == "same":
                struc_2d = {"disk": disk(6)}
                for i in range(np.shape(mask1_common)[2]):
                    mask1_common[:,:,i] = binary_dilation(mask1_common[:,:,i],struc_2d["disk"])
                seg1_common = (mask1_common) * seg1_common
                seg2_common = (mask1_common) * seg2_common
            elif mask == "different":
                seg1_common = (mask1_common) * seg1_common
                seg2_common = (mask2_common) * seg2_common
    elif component == "Cortical":
        mask1_common = cort1_common
        mask1_common[np.where(mask1_common != 0)] = 1
        mask2_common = cort2_common
        mask2_common[np.where(mask2_common != 0)] = 1
        if input_image == "Grayscale":
            if cort_surface == "both":
                if cort_mask == "same":
                    struc_2d = {"disk": disk(6)}
                    for i in range(np.shape(mask1_common)[2]):
                        mask1_common[:,:,i] = binary_dilation(mask1_common[:,:,i],struc_2d["disk"])
                    gray1_common = (mask1_common) * gray1_common
                    gray2_common = (mask1_common) * gray2_common
                    seg1_common = (mask1_common) * seg1_common
                    seg2_common = (mask1_common) * seg2_common
                elif cort_mask == "same_perio":
                    struc_2d = {"disk": disk(6)}
                    for i in range(np.shape(mask1_common)[2]):
                        mask1_common[:,:,i] = binary_dilation(mask1_common[:,:,i],struc_2d["disk"])
                    mask1_common = mask1_common-trab1_common
                    mask2_common = mask1_common-trab2_common
                    mask1_common[np.where(mask1_common == -1)] = 0
                    mask2_common[np.where(mask2_common == -1)] = 0
                    gray1_common = (mask1_common) * gray1_common
                    gray2_common = (mask2_common) * gray2_common
                    seg1_common = (mask1_common) * seg1_common
                    seg2_common = (mask2_common) * seg2_common
                elif cort_mask == "same_endo":
                    mask1_common = cort1_common + trab1_common
                    mask2_common = cort2_common + trab2_common
                    mask1_common = util.invert(flood_fill(mask1_common, (1, 1, 1), 1)) + 2 + mask1_common
                    mask2_common = util.invert(flood_fill(mask2_common, (1, 1, 1), 1)) + 2 + mask2_common
                    mask1_common[np.where(mask1_common != 0)] = 1    
                    mask2_common[np.where(mask2_common != 0)] = 1
                    struc_2d = {"disk": disk(6)}
                    for i in range(np.shape(mask1_common)[2]):
                        trab1_common_erod[:,:,i] = binary_erosion(trab1_common[:,:,i],struc_2d["disk"])
                    mask1_common = mask1_common-trab1_common
                    mask2_common = mask1_common-trab2_common
                    mask1_common[np.where(mask1_common == -1)] = 0
                    mask2_common[np.where(mask2_common == -1)] = 0
                    gray1_common = (mask1_common) * gray1_common
                    gray2_common = (mask2_common) * gray2_common
                    seg1_common = (mask1_common) * seg1_common
                    seg2_common = (mask2_common) * seg2_common
                elif cort_mask == "different":
                    gray1_common = (mask1_common) * gray1_common
                    gray2_common = (mask2_common) * gray2_common
                    seg1_common = (mask1_common) * seg1_common
                    seg2_common = (mask2_common) * seg2_common
            if cort_surface == "perio":
                mask1_common = cort1_common + trab1_common
                mask2_common = cort2_common + trab2_common
                mask1_common = util.invert(flood_fill(mask1_common, (1, 1, 1), 1)) + 2 + mask1_common
                mask2_common = util.invert(flood_fill(mask2_common, (1, 1, 1), 1)) + 2 + mask2_common
                mask1_common[np.where(mask1_common != 0)] = 1    
                mask2_common[np.where(mask2_common != 0)] = 1
                if cort_mask == "same_perio":
                    struc_2d = {"disk": disk(6)}
                    for i in range(np.shape(mask1_common)[2]):
                        mask1_common[:,:,i] = binary_dilation(mask1_common[:,:,i],struc_2d["disk"])
                        trab1_common[:,:,i] = binary_dilation(trab1_common[:,:,i],struc_2d["disk"])
                    mask1_common = mask1_common-trab1_common
                    mask1_common[np.where(mask1_common == -1)] = 0
                    gray1_common = (mask1_common) * gray1_common
                    gray2_common = (mask1_common) * gray2_common
                    seg1_common = (mask1_common) * seg1_common
                    seg2_common = (mask1_common) * seg2_common
                elif cort_mask == "different":
                    struc_2d = {"disk": disk(6)}
                    for i in range(np.shape(mask1_common)[2]):
                        trab1_common_dil[:,:,i] = binary_dilation(trab1_common[:,:,i],struc_2d["disk"])
                    mask1_common = mask1_common-trab1_common_dil
                    mask2_common = mask2_common-trab1_common_dil
                    gray1_common = (mask1_common) * gray1_common
                    gray2_common = (mask2_common) * gray2_common
                    seg1_common = (mask1_common) * seg1_common
                    seg2_common = (mask2_common) * seg2_common
            if cort_surface == "endo":
                mask1_common = cort1_common + trab1_common
                mask2_common = cort2_common + trab2_common
                mask1_common = util.invert(flood_fill(mask1_common, (1, 1, 1), 1)) + 2 + mask1_common
                mask2_common = util.invert(flood_fill(mask2_common, (1, 1, 1), 1)) + 2 + mask2_common
                mask1_common[np.where(mask1_common != 0)] = 1    
                mask2_common[np.where(mask2_common != 0)] = 1
                if cort_mask == "same_endo":
                    struc_2d = {"disk": disk(6)}
                    for i in range(np.shape(mask1_common)[2]):
                        mask1_common[:,:,i] = binary_erosion(mask1_common[:,:,i],struc_2d["disk"])
                        trab1_common[:,:,i] = binary_erosion(trab1_common[:,:,i],struc_2d["disk"])
                    mask1_common = mask1_common-trab1_common
                    gray1_common = (mask1_common) * gray1_common
                    gray2_common = (mask1_common) * gray2_common
                    seg1_common = (mask1_common) * seg1_common
                    seg2_common = (mask1_common) * seg2_common
                elif cort_mask == "different":
                    struc_2d = {"disk": disk(6)}
                    for i in range(np.shape(mask1_common)[2]):
                        mask1_common[:,:,i] = binary_erosion(mask1_common[:,:,i],struc_2d["disk"])
                    mask1_common = mask1_common-trab1_common
                    mask2_common = mask1_common-trab2_common
                    gray1_common = (mask1_common) * gray1_common
                    gray2_common = (mask2_common) * gray2_common
                    seg1_common = (mask1_common) * seg1_common
                    seg2_common = (mask2_common) * seg2_common
        elif input_image == "Binary":
            if cort_surface == "both":
                if cort_mask == "same":
                    struc_2d = {"disk": disk(6)}
                    for i in range(np.shape(mask1_common)[2]):
                        mask1_common[:,:,i] = binary_dilation(mask1_common[:,:,i],struc_2d["disk"])
                    seg1_common = (mask1_common) * seg1_common
                    seg2_common = (mask1_common) * seg2_common
                elif cort_mask == "same_perio":
                    struc_2d = {"disk": disk(6)}
                    for i in range(np.shape(mask1_common)[2]):
                        mask1_common[:,:,i] = binary_dilation(mask1_common[:,:,i],struc_2d["disk"])
                    mask1_common = mask1_common-trab1_common
                    mask2_common = mask1_common-trab2_common
                    mask1_common[np.where(mask1_common == -1)] = 0
                    mask2_common[np.where(mask2_common == -1)] = 0
                    seg1_common = (mask1_common) * seg1_common
                    seg2_common = (mask2_common) * seg2_common
                elif cort_mask == "same_endo":
                    mask1_common = cort1_common + trab1_common
                    mask2_common = cort2_common + trab2_common
                    mask1_common = util.invert(flood_fill(mask1_common, (1, 1, 1), 1)) + 2 + mask1_common
                    mask2_common = util.invert(flood_fill(mask2_common, (1, 1, 1), 1)) + 2 + mask2_common
                    mask1_common[np.where(mask1_common != 0)] = 1    
                    mask2_common[np.where(mask2_common != 0)] = 1
                    struc_2d = {"disk": disk(6)}
                    for i in range(np.shape(mask1_common)[2]):
                        trab1_common_erod[:,:,i] = binary_erosion(trab1_common[:,:,i],struc_2d["disk"])
                    mask1_common = mask1_common-trab1_common
                    mask2_common = mask1_common-trab2_common
                    mask1_common[np.where(mask1_common == -1)] = 0
                    mask2_common[np.where(mask2_common == -1)] = 0
                    seg1_common = (mask1_common) * seg1_common
                    seg2_common = (mask2_common) * seg2_common
                elif cort_mask == "different":
                    seg1_common = (mask1_common) * seg1_common
                    seg2_common = (mask2_common) * seg2_common
            if cort_surface == "perio":
                mask1_common = cort1_common + trab1_common
                mask2_common = cort2_common + trab2_common
                mask1_common = util.invert(flood_fill(mask1_common, (1, 1, 1), 1)) + 2 + mask1_common
                mask2_common = util.invert(flood_fill(mask2_common, (1, 1, 1), 1)) + 2 + mask2_common
                mask1_common[np.where(mask1_common != 0)] = 1    
                mask2_common[np.where(mask2_common != 0)] = 1
                if cort_mask == "same_perio":
                    struc_2d = {"disk": disk(6)}
                    for i in range(np.shape(mask1_common)[2]):
                        mask1_common[:,:,i] = binary_dilation(mask1_common[:,:,i],struc_2d["disk"])
                        trab1_common[:,:,i] = binary_dilation(trab1_common[:,:,i],struc_2d["disk"])
                    mask1_common = mask1_common-trab1_common
                    mask1_common[np.where(mask1_common == -1)] = 0
                    seg1_common = (mask1_common) * seg1_common
                    seg2_common = (mask1_common) * seg2_common
                elif cort_mask == "different":
                    struc_2d = {"disk": disk(6)}
                    for i in range(np.shape(mask1_common)[2]):
                        trab1_common_dil[:,:,i] = binary_dilation(trab1_common[:,:,i],struc_2d["disk"])
                    mask1_common = mask1_common-trab1_common_dil
                    mask2_common = mask2_common-trab1_common_dil
                    seg1_common = (mask1_common) * seg1_common
                    seg2_common = (mask2_common) * seg2_common
            if cort_surface == "endo":
                mask1_common = cort1_common + trab1_common
                mask2_common = cort2_common + trab2_common
                mask1_common = util.invert(flood_fill(mask1_common, (1, 1, 1), 1)) + 2 + mask1_common
                mask2_common = util.invert(flood_fill(mask2_common, (1, 1, 1), 1)) + 2 + mask2_common
                mask1_common[np.where(mask1_common != 0)] = 1    
                mask2_common[np.where(mask2_common != 0)] = 1
                if cort_mask == "same_endo":
                    struc_2d = {"disk": disk(6)}
                    for i in range(np.shape(mask1_common)[2]):
                        mask1_common[:,:,i] = binary_erosion(mask1_common[:,:,i],struc_2d["disk"])
                        trab1_common[:,:,i] = binary_erosion(trab1_common[:,:,i],struc_2d["disk"])
                    mask1_common = mask1_common-trab1_common
                    seg1_common = (mask1_common) * seg1_common
                    seg2_common = (mask1_common) * seg2_common
                elif cort_mask == "different":
                    struc_2d = {"disk": disk(6)}
                    for i in range(np.shape(mask1_common)[2]):
                        mask1_common[:,:,i] = binary_erosion(mask1_common[:,:,i],struc_2d["disk"])
                    mask1_common = mask1_common-trab1_common
                    mask2_common = mask1_common-trab2_common
                    seg1_common = (mask1_common) * seg1_common
                    seg2_common = (mask2_common) * seg2_common
    elif component == "Trabecular":
        mask1_common = trab1_common
        mask1_common[np.where(mask1_common != 0)] = 1
        mask2_common = trab2_common
        mask2_common[np.where(mask2_common != 0)] = 1
        if input_image == "Grayscale":
            if mask == "same":
                struc_2d = {"disk": disk(6)}
                for i in range(np.shape(mask1_common)[2]):
                    mask1_common[:,:,i] = binary_erosion(mask1_common[:,:,i],struc_2d["disk"])
                gray1_common = (mask1_common) * gray1_common
                gray2_common = (mask1_common) * gray2_common
                seg1_common = (mask1_common) * seg1_common
                seg2_common = (mask1_common) * seg2_common
            elif mask == "different":
                gray1_common = (mask1_common) * gray1_common
                gray2_common = (mask2_common) * gray2_common
                seg1_common = (mask1_common) * seg1_common
                seg2_common = (mask2_common) * seg2_common
        elif input_image == "Binary":
            if mask == "same":
                struc_2d = {"disk": disk(6)}
                for i in range(np.shape(mask1_common)[2]):
                    mask1_common[:,:,i] = binary_erosion(mask1_common[:,:,i],struc_2d["disk"])
                seg1_common = (mask1_common) * seg1_common
                seg2_common = (mask1_common) * seg2_common
            elif mask == "different":
                seg1_common = (mask1_common) * seg1_common
                seg2_common = (mask2_common) * seg2_common
    return gray1_common,gray2_common,seg1_common,seg2_common,mask1_common,mask2_common

def common_trim(seg1_common,seg2_common,cort1_common,cort2_common,trab1_common,trab2_common,histomorph):
    """
    Function for trimming the arrays to save space and time
    """
    # ===========================================================
    # Trimming the extra empty padding in all axes
    check_0 = ~(seg1_common==0).all((1,2))
    check_1 = ~(seg1_common==0).all((0,2))
    check_2 = ~(seg1_common==0).all((0,1))
    array_range_0 = np.where(check_0 == True)
    array_range_1 = np.where(check_1 == True)
    array_range_2 = np.where(check_2 == True)
    start0 = array_range_0[0][0]
    stop0 = array_range_0[0][-1]
    start1 = array_range_1[0][0]
    stop1 = array_range_1[0][-1]
    start2 = array_range_2[0][0]
    stop2 = array_range_2[0][-1]
    combined_img = histomorph[start0:stop0+1,start1:stop1+1,start2:stop2+1]
    seg1_common_tr = seg1_common[start0:stop0+1,start1:stop1+1,start2:stop2+1]
    seg2_common_tr = seg2_common[start0:stop0+1,start1:stop1+1,start2:stop2+1]
    cort1_common_tr = cort1_common[start0:stop0+1,start1:stop1+1,start2:stop2+1]
    cort2_common_tr = cort2_common[start0:stop0+1,start1:stop1+1,start2:stop2+1]
    trab1_common_tr = trab1_common[start0:stop0+1,start1:stop1+1,start2:stop2+1]
    trab2_common_tr = trab2_common[start0:stop0+1,start1:stop1+1,start2:stop2+1]

    return seg1_common_tr,seg2_common_tr, cort1_common_tr, cort2_common_tr, trab1_common_tr, trab2_common_tr,combined_img

def filter_remodeling(combined_img, seg1_common_tr, seg2_common_tr, cluster=5):
    """
    Function for removing small clusters (less than minimum cluster size) from bone remodeling labels (formation and resorption) 
    and generating a cleaned overlay map
    """
    def denoise(binary_img, min_size):
        labeled, _ = label(binary_img, background=0, connectivity=1, return_num=True)
        cleaned = remove_small_objects(labeled, min_size=min_size)
        return (cleaned > 0).astype(np.uint8)
    # Separate formation and resorption
    resorbed_bone = (combined_img == 500).astype(np.uint8)
    formed_bone = (combined_img == 300).astype(np.uint8)
    # Denoise both formation and resorption masks
    resorbed_denoised = denoise(resorbed_bone, cluster)
    formed_denoised = denoise(formed_bone, cluster)
    # Create overlay
    overlay = np.zeros_like(combined_img, dtype=np.uint16)
    overlay[seg1_common_tr == 1] = 600
    overlay[formed_denoised == 1] = 300
    overlay[resorbed_denoised == 1] = 500

    return overlay, formed_denoised, resorbed_denoised, combined_img, seg1_common_tr, seg2_common_tr

def binary_remodeling(combined_img,seg1_common_tr,seg2_common_tr):
    """
    Function for generating a cleaned overlay map for the binary method
    """
    resorbed_bone = np.zeros_like(combined_img)
    resorbed_bone = resorbed_bone.astype("uint8")
    formed_bone = np.zeros_like(resorbed_bone)
    constant_bone = np.zeros_like(resorbed_bone)
    # Creating separate arrays for resorption and formation
    resorbed_bone[np.where(combined_img == 500)] = 1
    formed_bone[np.where(combined_img == 300)] = 1

    overlay = np.zeros_like(formed_bone)
    overlay = overlay.astype(np.uint16)
    overlay[np.where(seg1_common_tr == 1)] = 600
    overlay[np.where(formed_bone == 1)] = 300
    overlay[np.where(resorbed_bone == 1)] = 500

    return overlay,formed_bone,resorbed_bone, combined_img, seg1_common_tr,seg2_common_tr


def calc_vol(formed,resorbed,base):
    """
    Function for calculating timelapse volume outcomes
    """
    MV = np.count_nonzero(formed)    # formed volume
    EV = np.count_nonzero(resorbed)    # resorbed volume
    BV = np.count_nonzero(base)  # baseline segmented volume
    MVBV = MV*100/BV
    EVBV = EV*100/BV
    
    return MV,EV,BV,MVBV,EVBV

def QC_tiff_export_vol(
    overlay, time1, time2, img_dir,
    patient_id=118, site="rad",
    scanner_generation=1, motion=0,
    registration="3D", input_image="Grayscale",
    erosion="no erosion", sigma=0.8,
    component="Total", mask="same",
    cort_mask="na", cort_surface="na",
    thresh=225, cluster=5,
):
    """
    Funciton for exporting Tiff stack containing timelapse outcome
    """
    label_mapping = {600: 255, 300: 172, 500: 117}
    overlay_mapped = np.copy(overlay)
    for label_value, gray_value in label_mapping.items():
        overlay_mapped[overlay == label_value] = gray_value

    overlay_mapped = overlay_mapped.astype(np.uint8)
    folder = os.path.join(
        img_dir,
        f"{component}_{cort_surface} bone"
    )
    os.makedirs(folder, exist_ok=True)

    tiff_name = (
        f"{patient_id}_{site}_"
        f"{time1}_{time2}_xct{scanner_generation}_{registration}_"
        f"{input_image}_{component}_{erosion}_{mask}_{cort_surface}_"
        f"{sigma}_{thresh}_{cluster}.tiff"
    )

    tiff_path = os.path.join(folder, tiff_name)
    imageio.mimwrite(tiff_path, overlay_mapped)

def flatten(seg1_common,seg2_common,cort1_common,cort2_common,trab1_common,trab2_common,histomorph):
    """
    Function for flattening the arrays. It finds slices to crop based on changes in cross-sectional area
    """
    mask1_common = cort1_common + trab1_common
    mask1_common = util.invert(flood_fill(mask1_common, (1, 1, 1), 1)) + 2 + mask1_common
    mask1_common[np.where(mask1_common != 0)] = 1    
    
    areas_list = np.array([0], int)
    # areas_list_orig = [0]
    for i in range(mask1_common.shape[2]):
        area = np.sum(mask1_common[:,:,i])
        areas_list = np.append(areas_list, area)
        areas_list = areas_list[areas_list != 0]   # Remove the empty slices
    areas_list_subt=np.ediff1d(areas_list, to_end=None, to_begin=0)      # Subtract consecutive areas
    res = np.divide(areas_list_subt, areas_list)   # Divide the difference in area by the area
    res_signs = np.copy(res)
    res_signs[res_signs < 0] = -1
    res_signs[res_signs > 0] = 1
    
    abs_res = abs(res)
    abs_med = np.quantile(abs_res, 0.5)   # Find the median change in area
    abs_div = abs_res/abs_med
    abs_div_bin = np.copy(abs_div)
    abs_div_bin[abs_div_bin < 5] = 0
    abs_div_bin[abs_div_bin >= 5] = 1    # Separate slices with changes in area 5 times larger than median
    signed_res_div = np.multiply(abs_div_bin,res_signs)
    # Identify large changes at the beginning and end
    begin_slices = [idx for idx, element in enumerate(signed_res_div) if element == 1]
    end_slices = [idx for idx, element in enumerate(signed_res_div) if element == -1]
    
    if len(begin_slices) != 0:
        begin_slice = begin_slices[-1] + 4
    else:
        begin_slice = 4
        
    if len(end_slices) != 0:
        end_slice = end_slices[0] - 4
    else:
        end_slice = int(mask1_common.shape[2])-4

    flat_region = str(begin_slice) + "_" + str(end_slice)
    
    seg1_common = seg1_common[:,:,begin_slice:end_slice]
    seg2_common = seg2_common[:,:,begin_slice:end_slice]
    cort1_common = cort1_common[:,:,begin_slice:end_slice]
    cort2_common = cort2_common[:,:,begin_slice:end_slice]
    trab1_common = trab1_common[:,:,begin_slice:end_slice]
    trab2_common = trab2_common[:,:,begin_slice:end_slice]
    histomorph = histomorph[:,:,begin_slice:end_slice]
    
    logger.info("Images are flattened to prevent any potential error!")
    
    return seg1_common,seg2_common,cort1_common,cort2_common,trab1_common,trab2_common,histomorph, flat_region, begin_slice,end_slice

def flatten_follow(
    seg1_common,seg2_common,cort1_common,cort2_common,
    trab1_common,trab2_common,histomorph,begin_slice,end_slice
):
    """
    Function for flattening the arrays using defined top and bottom slices
    """
    
    mask1_common = cort1_common + trab1_common
    mask1_common = util.invert(flood_fill(mask1_common, (1, 1, 1), 1)) + 2 + mask1_common
    mask1_common[np.where(mask1_common != 0)] = 1    

    flat_region = str(begin_slice) + "_" + str(end_slice)
    
    seg1_common = seg1_common[:,:,begin_slice:end_slice]
    seg2_common = seg2_common[:,:,begin_slice:end_slice]
    cort1_common = cort1_common[:,:,begin_slice:end_slice]
    cort2_common = cort2_common[:,:,begin_slice:end_slice]
    trab1_common = trab1_common[:,:,begin_slice:end_slice]
    trab2_common = trab2_common[:,:,begin_slice:end_slice]
    histomorph = histomorph[:,:,begin_slice:end_slice]
    
    logger.info("Images are flattened to prevent any potential error!")
    
    return seg1_common,seg2_common,cort1_common,cort2_common,trab1_common,trab2_common,histomorph, flat_region, begin_slice,end_slice

def run(
    gray1, gray2, seg1, seg2, cort1, cort2, trab1, trab2,
    header_gray1, header_gray2, header_seg1, header_seg2,
    header_cort1, header_cort2, header_trab1, header_trab2,
    time1, time2, csv_dir, image_repo_dir,
    Patient_ID=118, site="rad", scanner_generation=1,
    motion=0, registration="3D", input_image="Grayscale",
    erosion="no erosion", sigma=0.8, component="Total",
    mask="same", cort_mask="na", cort_surface="na",
    thresh=225, cluster=5, flatten_arrays=True, roi_begin=0, roi_ending=0
):
    """
    Function for running the timelapse calculations
    """
    logger.info(f"*********** Starting calculations for {component}_{cort_surface}_{sigma}_{thresh}_{cluster} ************")
    logger.info(f"++++++ Analyzing scans: {time1} and {time2} ++++++")
    
    # First, map all the image components to the common grid
    if input_image == "Grayscale":
        logger.info("Grayscale images will be used for calculations!")
    elif input_image == "Binary":
        logger.info("Binary images will be used for calculations!")
        
    gray1_common,gray2_common,seg1_common,seg2_common,cort1_common,cort2_common,trab1_common,trab2_common,header_gray1,header_gray2\
        = map2common_frame(gray1,gray2,seg1,seg2,cort1,cort2,trab1,trab2,header_gray1,header_gray2,header_seg1,header_seg2,\
                           header_cort1,header_cort2,header_trab1,header_trab2)
    # Apply Gaussian filter to the grayscale images
    if (sigma != 0) and (sigma != "na"): 
        gray1_common = gaussian_filter_3d(gray1_common,sigma)
        gray2_common = gaussian_filter_3d(gray2_common,sigma)
    
    # Remove the voxels at the edges of partial slices that are zero in one image due to aliasing
    if input_image == "Grayscale":
        gray1_common,gray2_common = gray_clean_up(gray1_common,gray2_common)

    # Mask the images to the region of interest
    gray1_common,gray2_common,seg1_common,seg2_common,mask1_common,mask2_common=\
        mask_images(gray1_common,gray2_common,seg1_common,seg2_common,cort1_common,cort2_common,trab1_common,trab2_common,input_image,component,mask,cort_mask,cort_surface)
    
    # Subtract the grayscale images to identify raw formation and resorption
    if input_image == "Grayscale":
        histomorph = gray_subtraction(gray1_common,gray2_common,header_gray1,header_gray2,thresh = thresh)
    elif input_image == "Binary":
        histomorph = binary_subtraction(seg1_common,seg2_common)

    # Trim the image components
    seg1_common,seg2_common,cort1_common,cort2_common,trab1_common,trab2_common,histomorph =\
        common_trim(seg1_common,seg2_common,cort1_common,cort2_common,trab1_common,trab2_common,histomorph)

    if flatten_arrays:
        if roi_begin == 0 and roi_ending == 0:
            # Flatten the images to prevent any potential error on the partial slices (only for total bone)
            seg1_common,seg2_common,cort1_common,cort2_common,trab1_common,trab2_common,histomorph, flat_region,begin_slice,end_slice = \
                flatten(seg1_common,seg2_common,cort1_common,cort2_common,trab1_common,trab2_common,histomorph)
        else: # If slices for flattening are given, apply them
            seg1_common,seg2_common,cort1_common,cort2_common,trab1_common,trab2_common,histomorph, flat_region,begin_slice,end_slice = \
                flatten_follow(seg1_common,seg2_common,cort1_common,cort2_common,trab1_common,trab2_common,histomorph,roi_begin,roi_ending)
    
    # Apply minimum cluster size (if minimum cluster size is zero, no changes will be made
    overlay,formed_denoised,resorbed_denoised, histomorph, seg1_common,seg2_common = filter_remodeling(histomorph,seg1_common,seg2_common,cluster)
    
    # Calculations
    logger.info("Doing calculations!")
    
    MV,EV,BV,MVBV,EVBV = calc_vol(formed_denoised,resorbed_denoised,seg1_common)
    
    QC_tiff_export_vol(overlay,time1,time2,image_repo_dir,Patient_ID,site,scanner_generation,motion,registration,input_image,erosion,sigma,component,mask,cort_mask,cort_surface,thresh, cluster)
  
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    
    # Prepare results for writing to CSV
    result_row = [
        Patient_ID, site, time1, time2, scanner_generation, motion, registration, input_image, erosion,
        sigma, component, mask, cort_mask, cort_surface, flat_region, thresh, cluster, MV, EV, BV, MVBV, EVBV, 
        date.today().isoformat(), current_time
    ]

    # Create a new CSV filename with date and time
    now = datetime.now()
    timestamp = now.strftime("%Y%m%d_%H%M%S")
    csv_filename = f"TL_results_{timestamp}.csv"
    csv_file_path = os.path.join(csv_dir, csv_filename)

    # Write the results to the new CSV file
    with open(csv_file_path, 'w', newline='') as file:
        writer = csv.writer(file)
        header = [
            "Patient_ID", "Site", "Time1", "Time2", "Scanner_Generation", "Motion", "Registration",
            "Input_Image", "Erosion", "Sigma", "Component", "Mask", "Cortical_Mask", "Cortical_Surface",
            "Flat_Region", "Threshold", "Cluster_Size", "MV", "EV", "BV", "MV/BV", "EV/BV", "Date", "Time"
        ]
        writer.writerow(header)

        # Then write the result row
        writer.writerow(result_row)
    
    logger.info("Time-lapse Computations Completed")

    return begin_slice, end_slice