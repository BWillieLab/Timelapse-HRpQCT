# Timelapse HR-pQCT

---

Python implementation of HR-pQCT timelapse analysis for quantifying bone formation and resorption:

***S. Hosseinitabatabaei et al, Non-invasive quantification of bone (re)modeling dynamics in adults with osteogenesis imperfecta treated with setrusumab using timelapse high-resolution peripheral-quantitative computed tomography, JBMR, 2025 \[(https://academic.oup.com/jbmr/article/40/3/348/7978263)\]***

In this implementation of timelapse HR-pQCT, noise reduction is achieved by Gaussian smoothing, instead of the common minimum cluster size filtering. Our results showed better performance (*i.e., smaller errors and noise, higher sensitivity*) compared to other methods.

The image below shows examples of the timelapse HR-pQCT using the Gaussian smoothing (*left column*) compared to minimum cluster size filtering (*details in the paper*):

![Visual examples of timelapse HR-pQCT with and without Gaussian smoothing](https://github.com/user-attachments/assets/64cea84f-8e35-4f7f-913f-9642596ce6b9)

This package can perform timelapse analysis using either grayscale or binary input images. However, the grayscale method is recommended, considering significantly lower errors. For the grayscale method, the following configurations are recommended (*details in the paper*):
- **Input image type:** Grayscale
- **Density threshold:** 200 mgHA/ccm
- **Minimum cluster size:** 0
- **Gaussian sigma:** 0.8 for XCT and 1.2 for XCT2

## Dependencies


- Python 3

- Python packages: `pip install packages listed below`
	- imageio==2.29.0
	- matplotlib==3.7.1
	- numpy==1.24.3
	- scipy==1.15.2
	- scikit-mage==0.20
	- vtk==9.2.6

	
## Dataset Preparation

This package requires the following Scanco-generated AIM files of scan pairs as inputs, aligned and segmented using IPL.

The scan pairs that need to be analyzed should be stored in a parent folder, each in a subfolder that preferably indicates anatomical site in the subfolder name.

```
parent_folder/
├── ID0001_radius
│   	├── 1_gray1.aim
│ 	├── 2_gray2.aim
│ 	├── 3_cortmask1.aim
│ 	├── 4_cortmask2.aim
│ 	├── 5_trabmask1.aim
│ 	├── 6_trabmask2.aim
│ 	├── 7_seg1.aim
│ 	└── 8_seg2.aim
└── ID0001_tibia
    	├── 1_gray1.aim
 	├── 2_gray2.aim
	├── 3_cortmask1.aim
	├── 4_cortmask2.aim
 	├── 5_trabmask1.aim
 	├── 6_trabmask2.aim
 	├── 7_seg1.aim
	└── 8_seg2.aim
```


## Installation & Usage

### 🔧 Installation

To install this package manually:

1. Download or clone the repository:
   ```bash
   git clone https://github.com/BWillieLab/Timelapse_HRpQCT.git
   ```

2. Go to the project folder:
   ```bash
   cd Timelapse_HRpQCT
   ```

3. Run your Python script as usual:
   ```bash
   python hrpqct_tl_run.py -c options/tl_config.json
   ```

> 💡 Make sure you have the required Python libraries (listed at the top of the document) are installed. If something is missing, install it using:
```bash
pip install library-name
```



### ▶️ Running the Analyses

**STEP 1:** Open the `tl_config.json` file using a text editor. You need to specify paths for `parent_folder`, `image_repo`, and `csv_dir`. Other parameters have default values (described below):

```json
{
  "parent_folder": "/path/to/parent_folder",
  "parameters": {
    "timepoint1": "Bsl",
    "timepoint2": "6M",
    "sigma": 1.2,
    "component": "Total",
    "mask": "same",
    "cort_mask": "na",
    "cort_surface": "na",
    "thresh": 200,
    "cluster": 0,
    "roi_begin": 0,
    "roi_ending": 0,
    "image_repo": "/path/to/image_repo",
    "csv_dir": "/path/to/csv_file"
  }
}
```

### 📝 Parameter Descriptions

| Key | Description |
|-----|-------------|
| `parent_folder` | Folder containing subfolders with scan pairs |
| `timepoint1` | Name of the first timepoint (e.g., Baseline) |
| `timepoint2` | Name of the second timepoint (e.g., 6-month follow-up) |
| `sigma` | Gaussian filter sigma (use 0.8 for XCT, 1.2 for XCT2) |
| `component` | Bone compartment to analyze (`Total`, `Cortical`, `Trabecular`) |
| `mask` | Whether the periosteal mask is the same (`same`) or different |
| `cort_mask` | Same or different cortical masks (`same`, `different`, or `na` if not used) |
| `cort_surface` | Cortical surface type (`perio`, `endo`, `both`, or `na`) |
| `thresh` | Density threshold for segmentation (e.g., 200) |
| `cluster` | Minimum cluster size to retain (e.g., 0 = no filtering) |
| `roi_begin` | Number of top slices to crop |
| `roi_ending` | Number of bottom slices to crop |
| `image_repo` | Directory to save TIFF stacks for QC and visualization |
| `csv_dir` | Directory to save the resulting CSV files |

**STEP 2:** Make sure you are inside the `code/` folder of the project.

**STEP 3:** Run the analysis with:
```bash
python hrpqct_tl_run.py -c options/tl_config.json
```


## Contact 
For any questions, contact <mahdi.tabatabaei@mail.mcgill.ca>
