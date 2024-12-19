# Lung Segmentation on CT Scans

## Project Overview

This project focuses on segmenting lung structures from CT scans, with both **automatic** and **manual** approaches. The workflow is designed to evaluate and compare segmentation techniques, calculate lung areas and volumes, and analyze segmentation performance using metrics like the **Dice coefficient**, **Average Surface Distance (ASD)**, and **local curvature differences**. 

### Objectives

1. **Segment the lungs** from CT scans and calculate the cross-sectional areas for axial and coronal slices.
2. **Quantify total lung volume** and height based on segmented regions.
3. **Introduce noise** to evaluate segmentation robustness.
4. **Compare manual and automatic segmentations** using performance metrics.
5. Provide **3D visualization** of segmented lung structures.

## Features

- **Interactive slice selection** from CT volumes for segmentation.
- **Region of Interest (ROI) cropping** to optimize processing.
- **Automatic segmentation** using Otsu's thresholding and morphological operations.
- **Manual segmentation** with user-defined contours.
- **Comparison metrics**:
  - Dice coefficient.
  - Average Surface Distance (ASD).
  - Relative area difference.
  - Local curvature analysis.
- **Noise evaluation**: Add Gaussian or pepper noise to assess segmentation robustness.
- **Volume estimation** and 3D visualization.

---

## Workflow

### 1. Data Preparation

- Input: 3D CT volume in DICOM format.
- A **user interface** enables browsing through slices and selecting regions for analysis.

### 2. Segmentation

#### Automatic Segmentation
- Converts the image to grayscale.
- Applies Otsu's method for automatic thresholding.
- Uses morphological operations to refine the segmented regions.

#### Manual Segmentation
- Users select lung contours by clicking points.
- Masks are created based on the selected points.

### 3. Comparison

Both segmentation techniques are compared using:
- Dice coefficient.
- Average Surface Distance (ASD).
- Relative area difference.
- Curvature analysis.

### 4. Noise Evaluation

- Add Gaussian or pepper noise to the dataset.
- Re-run segmentation workflows to analyze the robustness.

### 5. Visualization

- Overlay automatic and manual contours for comparison.
- Plot curvature differences.
- 3D visualization of segmented volumes using `volumeViewer`.

---

## Usage

### Prerequisites

- MATLAB with the following toolboxes:
  - Image Processing Toolbox.
  - Statistics and Machine Learning Toolbox (optional).
- A DICOM dataset. [Example dataset from Kitware](https://data.kitware.com/#collection/579787098d777f1268277a27/folder/5a9dc8f78d777f06857860fd).

### Steps to Run

1. Clone this repository:
   ```bash
   git clone <repository-url>
   cd <repository-directory>
   ```
2. Ensure the dataset path is updated in the scripts.
3. Run the desired segmentation script in MATLAB:
   ```matlab
   ctLungSegmentationAxial;
   ```
   or
   ```matlab
   ctLungSegmentationCoronal;
   ```

### Functions Overview

| Function                     | Description                                                                 |
|------------------------------|-----------------------------------------------------------------------------|
| `readDCMfolder`              | Reads DICOM files and generates a 3D volume.                               |
| `visualizeAndSelectSlice`    | Enables interactive selection of slices.                                   |
| `segmentLungAutomatically`   | Performs automatic lung segmentation.                                      |
| `segmentLungManually`        | Enables manual lung segmentation by user-defined contours.                 |
| `calculateDice`              | Computes the Dice coefficient.                                             |
| `computeAverageSurfaceDist`  | Calculates Average Surface Distance between contours.                      |
| `computeLocalCurvature`      | Computes local curvatures for contour comparison.                          |
| `overlayContours`            | Overlays manual and automatic contours for visual inspection.              |
| `applySegmentationToAllSlices`| Segments and processes all slices in a volume, calculates areas and volume.|
| `readDCMfolder`              | Loads DICOM files, constructs the 3D volume, and extracts metadata.        |

---

## Results

### Metrics

| Metric                         | Axial Planes   | Coronal Planes |
|--------------------------------|----------------|----------------|
| **Largest Cross-sectional Area** | ~20,560 mm²    | ~19,500 mm²    |
| **Total Lung Volume**           | ~2.41 L        | ~2.02 L        |
| **Height**                      | ~16.5 cm       | ~15.75 cm      |

### Observations

- **Manual vs Automatic Comparison**:
  - Dice coefficient: ~0.95.
  - ASD: ~3-4 mm (axial), ~1-2 mm (coronal).
  - Relative area difference: ~6-7%.
  - Local curvature differences: Minimal.
- **Noise Impact**:
  - Increased noise degrades segmentation accuracy.
  - Pepper noise has a more pronounced effect on edge detection.

---

## Future Work

- Explore **machine learning** or **deep learning models** for segmentation (e.g., U-Net).
- Apply the workflow to diseased lungs for further evaluation.
- Extend comparison to additional metrics or larger datasets.

---

## References

1. Candemir, S., & Antani, S. (2019). A review on lung boundary detection in chest X-rays. *International Journal of Computer Assisted Radiology and Surgery*, 14(4), 563–576.
2. Mansoor, A., et al. (2015). Segmentation and Image Analysis of Abnormal Lungs at CT: Current Approaches, Challenges, and Future Trends. *RadioGraphics*, 35(4), 1056–1076.

---

## Authors

- [Maë Lebedel](https://github.com/maelebedel)

Feel free to contribute to the project by submitting issues or pull requests!
