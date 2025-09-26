"""
Script Name: hrpqct_tl_call.py
Description: This script:
                - Organizes the AIM files in the folder of each scan pair
                - Passes the scan pair files and configs into the hrpqct_tl_comp.py script for processing
"""

from utils.logger_setup import logger
import os
from pathlib import Path
from timelapse_analysis import hrpqct_tl_comp as tl
from utils import aim_reader_v1 as aim


def timelapse_call(dirname, config_params):

    timepoint1 = config_params.get("timepoint1", "Bsl")
    timepoint2 = config_params.get("timepoint2", "6M")
    timepoint3 = config_params.get("timepoint3", "12M")
    timepoint4 = config_params.get("timepoint4", "18M")
    timepoint5 = config_params.get("timepoint5", "24M")
    sigma = config_params.get("sigma", 1.2)
    component = config_params.get("component", "Total")
    mask = config_params.get("mask", "same")
    cort_mask = config_params.get("cort_mask", "na")
    cort_surface = config_params.get("cort_surface", "na")
    thresh = config_params.get("thresh", 200)
    cluster = config_params.get("cluster", 0)
    roi_beg = config_params.get("roi_begin", 0)
    roi_end = config_params.get("roi_ending", 0)
    image_repo_dir = config_params.get("image_repo")
    csv_dir = config_params.get("csv_dir")
    
    all_files = os.listdir(dirname)

    input_files_norder = [f for f in all_files if "AIM" in f]

    input_files = []
    for fl in range(1, 33):
        for fll in input_files_norder:
            under = fll.index("_")
            if fll[:under] == str(fl):
                input_files.append(fll)

    logger.info(input_files)

    dirnamep = Path(dirname)
    
    # Create full path using Path object for cross-platform compatibility
    input_files = [str(dirnamep / f) for f in input_files]

    # Sorting scans by date, and finding which timepoints are missing
    # First, read the grayscale images by number to extract their scan dates
    init_gray_order = [i for i in range(8)]
    init_cort_order = [i + 8 for i in range(8)]
    init_trab_order = [i + 16 for i in range(8)]
    init_seg_order =  [i + 24 for i in range(8)]
    
    aim_gray1 = input_files[0]
        # Load images and headers...
    gray1, header_gray1 = aim.aim2np(aim_gray1)
    Pat_ID = header_gray1["Patient_ID"]
    Site = "Radius" if "rad" in dirname.lower() else "Tibia"
    scanner = 1 if header_gray1["XCT_gen"] == 1 else 2

    logger.info(f"The scans are taken from this folder: {dirname}")

    # # Now, findings and sorting the input files
    aim_gray1, aim_gray2, aim_gray3, aim_gray4, aim_gray5, aim_gray6, aim_gray7, aim_gray8 = [input_files[i] for i in init_gray_order]
    aim_cort1, aim_cort2, aim_cort3, aim_cort4, aim_cort5, aim_cort6, aim_cort7, aim_cort8 = [input_files[i] for i in init_cort_order]
    aim_trab1, aim_trab2, aim_trab3, aim_trab4, aim_trab5, aim_trab6, aim_trab7, aim_trab8 = [input_files[i] for i in init_trab_order]
    aim_seg1, aim_seg2, aim_seg3, aim_seg4, aim_seg5, aim_seg6, aim_seg7, aim_seg8         = [input_files[i] for i in init_seg_order]

    logger.info("reading AIM files")
    gray2,header_gray2 = aim.aim2np(aim_gray2)
    gray3,header_gray3 = aim.aim2np(aim_gray3)
    gray4,header_gray4 = aim.aim2np(aim_gray4)
    gray5,header_gray5 = aim.aim2np(aim_gray5)
    gray6,header_gray6 = aim.aim2np(aim_gray6)
    gray7,header_gray7 = aim.aim2np(aim_gray7)
    gray8,header_gray8 = aim.aim2np(aim_gray8)

    cort1,header_cort1 = aim.aim2np(aim_cort1)
    cort2,header_cort2 = aim.aim2np(aim_cort2)
    cort3,header_cort3 = aim.aim2np(aim_cort3)
    cort4,header_cort4 = aim.aim2np(aim_cort4)
    cort5,header_cort5 = aim.aim2np(aim_cort5)
    cort6,header_cort6 = aim.aim2np(aim_cort6)
    cort7,header_cort7 = aim.aim2np(aim_cort7)
    cort8,header_cort8 = aim.aim2np(aim_cort8)

    trab1,header_trab1 = aim.aim2np(aim_trab1)
    trab2,header_trab2 = aim.aim2np(aim_trab2)
    trab3,header_trab3 = aim.aim2np(aim_trab3)
    trab4,header_trab4 = aim.aim2np(aim_trab4)
    trab5,header_trab5 = aim.aim2np(aim_trab5)
    trab6,header_trab6 = aim.aim2np(aim_trab6)
    trab7,header_trab7 = aim.aim2np(aim_trab7)
    trab8,header_trab8 = aim.aim2np(aim_trab8)
    
    seg1,header_seg1 = aim.aim2np(aim_seg1)
    seg2,header_seg2 = aim.aim2np(aim_seg2)
    seg3,header_seg3 = aim.aim2np(aim_seg3)
    seg4,header_seg4 = aim.aim2np(aim_seg4)
    seg5,header_seg5 = aim.aim2np(aim_seg5)
    seg6,header_seg6 = aim.aim2np(aim_seg6)
    seg7,header_seg7 = aim.aim2np(aim_seg7)
    seg8,header_seg8 = aim.aim2np(aim_seg8)

    logger.info("Performing the computations (Bsl and 6M)")
    tl.run(
        gray1, gray2, seg1, seg2, cort1, cort2, trab1, trab2,
        header_gray1, header_gray2, header_seg1, header_seg2,
        header_cort1, header_cort2, header_trab1, header_trab2,
        timepoint1, timepoint2,csv_dir, image_repo_dir,
        Patient_ID=Pat_ID, site=Site, scanner_generation=scanner,
        sigma=sigma, component=component, mask=mask,
        cort_mask=cort_mask, cort_surface=cort_surface,
        thresh=thresh, cluster=cluster,
        flatten_arrays=True, roi_begin=roi_beg, roi_ending=roi_end
    )

    logger.info("Performing the computations (6M and 12M)")
    tl.run(
        gray3, gray4, seg3, seg4, cort3, cort4, trab3, trab4,
        header_gray3, header_gray4, header_seg3, header_seg4,
        header_cort3, header_cort4, header_trab3, header_trab4,
        timepoint2, timepoint3, csv_dir, image_repo_dir,
        Patient_ID=Pat_ID, site=Site, scanner_generation=scanner,
        sigma=sigma, component=component, mask=mask,
        cort_mask=cort_mask, cort_surface=cort_surface,
        thresh=thresh, cluster=cluster,
        flatten_arrays=True, roi_begin=roi_beg, roi_ending=roi_end
    )
    
    logger.info("Performing the computations (12M and 18M)")
    tl.run(
        gray5, gray6, seg5, seg6, cort5, cort6, trab5, trab6,
        header_gray5, header_gray6, header_seg5, header_seg6,
        header_cort5, header_cort6, header_trab5, header_trab6,
        timepoint3, timepoint4, csv_dir, image_repo_dir,
        Patient_ID=Pat_ID, site=Site, scanner_generation=scanner,
        sigma=sigma, component=component, mask=mask,
        cort_mask=cort_mask, cort_surface=cort_surface,
        thresh=thresh, cluster=cluster,
        flatten_arrays=True, roi_begin=roi_beg, roi_ending=roi_end
    )

    logger.info("Performing the computations (18M and 24M)")
    tl.run(
        gray7, gray8, seg7, seg8, cort7, cort8, trab7, trab8,
        header_gray7, header_gray8, header_seg7, header_seg8,
        header_cort7, header_cort8, header_trab7, header_trab8,
        timepoint4, timepoint5, csv_dir, image_repo_dir,
        Patient_ID=Pat_ID, site=Site, scanner_generation=scanner,
        sigma=sigma, component=component, mask=mask,
        cort_mask=cort_mask, cort_surface=cort_surface,
        thresh=thresh, cluster=cluster,
        flatten_arrays=True, roi_begin=roi_beg, roi_ending=roi_end
    )
