# STIX-Solar-Orbiter

# STIX Spectrum Fitting Tool in Python

This project provides tools to fit spectra measured by the STIX instrument (Solar Orbiter) using custom spectral models and the instrument response matrix (SRM). It features an interactive Tkinter-based interface for loading, modeling, and comparing data.

## 📦 Features

- Load STIX FITS files.
- Choose from various spectral models: PowerLaw1D, BrokenPowerLaw1D, V_TH, etc.
- Forward folding using the SRM matrix.
- Automatic model fitting.
- Interactive visualization of results (flux, rate, counts).
- Support for statistical error propagation.

## 🖥️ Interface

The Tkinter GUI allows you to:
- Load FITS files (spectra and SRM).
- Select energy intervals.
- Add and configure spectral models.
- Perform fitting and display results.


## ⚙️ Installation

### Requirements

- Python 3.9+
- Astropy
- Numpy
- Scipy
- Matplotlib
- Tkinter

### Setup

    git clone https://github.com/Assamoi21/STIX-Solar-Orbiter.git
    cd STIX-Solar-Orbiter
    pip install -r requirements.txt


## 📁 Data

Example FITS files (spectra and SRM) are included in the repository. You can also download them from the official Solar Orbiter data sources.

## 📜 Licence

This project is licensed under the MIT License.

## 👨‍🔬 Authors

    Abdallah Hamini, Assamoua Koman

    Contact : abdallah.hamini@obspm.fr


