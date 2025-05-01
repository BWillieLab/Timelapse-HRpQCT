import logging
import os
from datetime import datetime

logger = logging.getLogger("TimeLapseLogger")

def setup_logger(parent_folder=".", log_prefix="time_lapse_log"):
    if logger.handlers:
        return logger
    os.makedirs(parent_folder, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_filename = f"{log_prefix}_{timestamp}.txt"
    log_file_path = os.path.join(parent_folder, log_filename)
    logger.setLevel(logging.DEBUG)

    file_handler = logging.FileHandler(log_file_path, mode='w')
    file_handler.setLevel(logging.DEBUG)

    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    file_handler.setFormatter(formatter)

    logger.addHandler(file_handler)

    return logger
