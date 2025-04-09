Phase-Amplitude Coupling (PAC) Analysis Toolkit
This repository contains MATLAB code for calculating Phase-Amplitude Coupling (PAC) metrics as described in our STAR Protocols paper and previously implemented in Aliramezani et al., 2024. 

📦 What's Included
Core Functions:
•	Main_PAC.m: The main script that loads data and calculates PAC measures.
•	PACcalculator_MI.m: Computes PAC using modulation index (MI) method (Tort et al., 2010).
•	Morlet_Wavelet.m: Performs wavelet decomposition for phase/amplitude extraction
•	Notch_Filter.m: Removes power line noise (50/60 Hz).
•	Remove_Artifact.m: Identifies artifact-contaminated trials.

Sample Data
Three sample datasets from our associated paper (Aliramezani et al., 2024) are included to demonstrate functionality:
1.	SampleData1.mat.
2.	SampleData2.mat.
3.	SampleData3.mat.

🚀 Getting Started
Basic Usage:
1.	Clone this repository or download the files.
2.	Open MATLAB and navigate to the repository folder.
3.	Run Main_PAC.m to:
o	Load sample data.
o	Calculate PAC metrics.
o	Generate results.
Customizing for Your Data. To analyze your own data:
1.	Replace the sample data file name in Main_PAC.m.
2.	Adjust parameters in the "User-defined parameters" section.

🔍 Method Overview
The code implements:
1.	Artifact rejection (threshold-based).
2.	Notch filtering (50/60 Hz line noise removal).
3.	Wavelet decomposition (Morlet wavelets).
4.	PAC calculation using Modulation Index (MI).
For full methodological details, please refer to our STAR Protocols paper.

🤝 How to Cite
If you use this code in your research, please cite both:
1.	Our STAR Protocols paper (citation available upon publication):
Aliramezani, M., Farrokhi, A., Constantinidis, C., and Daliri, M.R. (2025). Protocol for phase-amplitude coupling analysis in local field potentials to investigate neural oscillation dynamics.
2.	Our associated paper:
Aliramezani, M., Farrokhi, A., Constantinidis, C., and Daliri, M.R. (2024). Delta-alpha/beta coupling as a signature of visual working memory in the prefrontal cortex. iScience 27, 110453. https://doi.org/10.1016/j.isci.2024.110453.
😉Dream biG ...
