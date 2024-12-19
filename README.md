# Lung Segmentation on CT Scans

## Project Overview

This project focuses on segmenting lung structures from CT scans, with both **automatic** and **manual** approaches. The workflow is designed to evaluate and compare segmentation techniques, calculate lung areas and volumes, and analyze segmentation performance using metrics like the **Dice coefficient**, **Average Surface Distance (ASD)**, **Relative Area Difference**, and **Local Curvature Analysis**. 

This project also aims to highlight the challenges associated with accurately segmenting lung boundaries, which often involve complex and variable edge structures.

### Objectives

1. **Segment the lungs** from CT scans and calculate the cross-sectional areas for axial and coronal slices.
2. **Quantify total lung volume** and height based on segmented regions.
3. **Introduce noise** to evaluate segmentation robustness.
4. **Compare manual and automatic segmentations** using performance metrics.
5. Provide **3D visualization** of segmented lung structures.
6. **Highlight the challenges** of segmenting lung boundaries.

## Features

- **Interactive slice selection** from CT volumes for segmentation.
- **Region of Interest (ROI) cropping** to optimize processing.
- **Automatic segmentation** using Otsu's thresholding and morphological operations.
- **Manual segmentation** with user-defined contours.
- **Comparison metrics**:
  - Dice coefficient.
  - Average Surface Distance (ASD).
  - Relative Area Difference.
  - Local Curvature Analysis.
- **Noise evaluation**: Add Gaussian or pepper noise to assess segmentation robustness.
- **Volume estimation** and 3D visualization.

---

## Workflow

### 1. Data Preparation

- Input: 3D CT volume in DICOM format.
- A **user interface** enables browsing through slices and selecting regions for analysis.
- The dataset contains partial lung images; the scans do not fully capture the entire lung structures. This limitation explains the relatively low lung volume and area measurements compared to typical full-lung values.
- **Dataset Details**:
  - Source: [Kitware DICOM Dataset](https://data.kitware.com/#collection/579787098d777f1268277a27/folder/5a9dc8f78d777f06857860fd).
  - Steps to access:
    1. Register on the Kitware Data platform.
    2. Navigate to `patient0/T_0/CT` and download the DICOM files.
    3. Use the provided `readDCMfolder` function to load the dataset into a MATLAB 3D volume stack.

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

#### Metrics Explanation:
- **Dice coefficient**:
  - Measures the overlap between two segmentations.
  - **1**: Perfect overlap.
  - **0**: No overlap.
- **Average Surface Distance (ASD)**:
  - Measures the average distance between the contours of two segmentations.
- **Relative Area Difference**:
  - Compares the areas of two segmentations as a percentage.
  - Calculated as \( \frac{|A_{manual} - A_{auto}|}{A_{manual}} \times 100 \% \).
- **Local Curvature Analysis**:
  - Compares contour curvatures to identify differences in smoothness or shape.
  - Useful for detecting regions with higher curvature (e.g., sharp edges).

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

- **MATLAB** with the following toolbox:
  - Image Processing Toolbox (required).
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

Each part of the project contains distinct functions, but the following are commonly used:

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

---

## Results

### Metrics

| Metric                         | Axial Planes                     | Coronal Planes                   |
|--------------------------------|-----------------------------------|-----------------------------------|
| **Largest Cross-sectional Area** | ~20,560 mm² (Axial only)         | N/A                              |
| **Total Lung Volume**           | ~2.41 L                          | ~2.02 L                          |
| **Height**                      | ~16.5 cm                         | ~15.75 cm                        |
| **Dice coefficient**            | ~0.95                            | ~0.93                            |
| **ASD**                         | ~3-4 mm                          | ~1-2 mm                          |
| **Relative Area Difference**    | ~6-7%                            | ~5-6%                            |
| **Local Curvature Differences** | Minimal                          | Minimal                          |

### Observations

- **Manual vs Automatic Comparison**:
  - Axial planes tend to show slightly higher Dice scores and ASD values due to the alignment and partial visualization of structures.
  - Coronal planes offer smoother contours and lower ASD due to better continuity along the scan.
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
