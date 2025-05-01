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
    for fl in range(1, 9):
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
    init_gray_order = [i for i in range(2)]
    init_cort_order = [i + 2 for i in range(2)]
    init_trab_order = [i + 4 for i in range(2)]
    init_seg_order = [i + 6 for i in range(2)]
    
    aim_gray1 = input_files[0]
        # Load images and headers...
    gray1, header_gray1 = aim.aim2np(aim_gray1)
    Pat_ID = header_gray1["Patient_ID"]
    Site = "Radius" if "rad" in dirname.lower() else "Tibia"
    scanner = 1 if header_gray1["XCT_gen"] == 1 else 2

    logger.info(f"The scans are taken from this folder: {dirname}")

    # # Now, findings and sorting the input files
    aim_gray1, aim_gray2 = [input_files[i] for i in init_gray_order]
    aim_cort1, aim_cort2 = [input_files[i] for i in init_cort_order]
    aim_trab1, aim_trab2 = [input_files[i] for i in init_trab_order]
    aim_seg1, aim_seg2 = [input_files[i] for i in init_seg_order]

    logger.info("reading AIM files")
    gray2,header_gray2 = aim.aim2np(aim_gray2)
    cort1,header_cort1 = aim.aim2np(aim_cort1)
    cort2,header_cort2 = aim.aim2np(aim_cort2)
    trab1,header_trab1 = aim.aim2np(aim_trab1)
    trab2,header_trab2 = aim.aim2np(aim_trab2)
    seg1,header_seg1 = aim.aim2np(aim_seg1)
    seg2,header_seg2 = aim.aim2np(aim_seg2)

    logger.info("Performing the computations")
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