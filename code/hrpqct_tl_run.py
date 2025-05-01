"""
Script Name: hrpqct_tl_run.py
Description: This script:
                - Collects the configurations
                - Creates the logger to log the analyses
                - Identify the scan pairs in the parent folder to be processed
                - Passes the scan pairs into the hrpqct_tl_call.py script for preparation
"""

import json
import os
import argparse
from utils import logger_setup
from timelapse_analysis import hrpqct_tl_call as tlc

def main(config_path):
    # Load json file containing the paramters
    with open(config_path, 'r') as f:
        config = json.load(f)

    parent_folder = config["parent_folder"]
    params = config["parameters"]
    logger = logger_setup.setup_logger(parent_folder)

    logger.info("STARTING THE BATCH TIMELAPSE ANALYSES")

    # Log config options here
    logger.info(f"Loaded configuration from {config_path}")
    logger.debug(f"Configuration content: {json.dumps(config, indent=2)}")

    # Going through subfolders
    for root, dirs, files in os.walk(parent_folder):
        for d in dirs:
            full_path = os.path.join(root, d)
            if any("AIM" in f for f in os.listdir(full_path)):
                logger.info(f"Processing folder: {full_path}")
                tlc.timelapse_call(full_path, params)

    logger.info("TIMELAPSE ANALYSES COMPLETED ON ALL SCAN SETS")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run batch timelapse analysis with a config file.")
    parser.add_argument(
        "-c", "--config",
        type=str,
        required=True,
        help="Path to the JSON config file."
    )
    args = parser.parse_args()
    main(args.config)
