## Overview

This repository contains the code associated with the paper **"On-the-fly spectral unmixing based on Kalman filtering"**. In this paper, we propose performing spectral unmixing on a spectrum-by-spectrum basis.

## Contents

### Data Test Folder

- **Pure Spectra:**

   - `pureSpectra_matrix.mat`: Pure spectra used in the paper.

- **Concentration Matrices:**

  - `concentration_matrix_SD1_P1.mat`: A concentration matrix generated as described in data sets SD1, protocol P1.
  - `concentration_matrix_SD1_P2.mat`: The associated "essential concentration" matrix.
  - `concentration_matrix_SD2_P1.mat`: A concentration matrix generated as described in data sets SD2, protocol P1.
  - `concentration_matrix_SD2_P2.mat`: The associated "essential concentration" matrix.
  - `concentration_matrix_SubspaceTracking.mat`: A concentration matrix for evaluating the SubspaceTracking algorithm, where 3 components are present from time index 1 to 50, an additional one appears from time index 51 to 1000, and a final one appears from time index 1001 to 4000.

### KFOSU Folder
   - Contains KF-OSU code.

### VCAandSUNSAL Folder
   - Includes the Vertex Component Analysis (VCA) and Sparse Unmixing by Variable Splitting and Augmented Lagrangian (SUnSAL) code.

## Usage

To run the code, execute `Test.m`. This script demonstrates the use of KF-OSU with an example that uses one of the concentration matrices listed above. You can adjust the signal-to-noise ratio as well as other parameters, as explained in `Test.m` and `KFOSU.m`.

## Citation

If you use this code, please cite the paper as follows:

```bibtex
@article{kouakou2024fly,
  title={On-the-fly spectral unmixing based on Kalman filtering},
  author={Kouakou, Hugues and de Morais Goulart, Jos{\'e} Henrique and Vitale, Raffaele and Oberlin, Thomas and Rousseau, David and Ruckebusch, Cyril and Dobigeon, Nicolas},
  journal={Chemometrics and Intelligent Laboratory Systems},
  pages={105252},
  year={2024},
  publisher={Elsevier}
}
