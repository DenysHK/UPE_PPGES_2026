[![DOI](https://zenodo.org/badge/1184710780.svg)](https://doi.org/10.5281/zenodo.19204383)

# Modeling and Simulation of Dry Reforming of Biogas and Fischer–Tropsch Synthesis

This repository contains the computational models and scripts developed within a Master's dissertation in Systems Engineering (PPGES – UPE).

## 🔬 Scope

The repository includes:

* Dry Reforming of Biogas (RSB) modeling
* Syngas production analysis
* Fischer–Tropsch synthesis (SFT) modeling
* Product distribution via Anderson–Schulz–Flory (ASF)
* Climate data analysis (C3S)

## 📁 Repository Structure

### `dwsim/`

Base simulation models developed in DWSIM:

* `RSB_Benguerba2015_BaseModel.dwxmz`
* `SFT_Pandey2021_BaseModel.dwxmz`

### `python/`

Supporting scripts:

* `ASF_Alpha_Distribution_Plot.py`
* `C3S_GlobalTemperature_AnnualAbove1p5C_Plot.py`
* `SFT_Pandey2021_LHHW_MassBalance_DWSIM.py`

## ⚙️ Requirements

* Python 3.x
* DWSIM
* Standard scientific libraries (NumPy, Matplotlib, Pandas)

## 📌 Reproducibility

All models and scripts are provided to ensure reproducibility of the results presented in the dissertation.

## 📚 Citation

If you use this repository, please cite:

KABBAZ, Denys Haluli. Modeling and simulation of dry reforming of biogas and Fischer–Tropsch synthesis. 2026. Software. DOI: (to be inserted)

## 📄 License

To be defined (recommended: MIT License)
