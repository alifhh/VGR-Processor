# VGR Processor (v1.0.0)

**VGR Processor** is a standalone, zero-dependency Python-based software designed for processing raw multicomponent volcanic gas analyzer (Multi-Gas) data. It provides an intuitive graphical user interface (GUI) to handle complex signal processing workflows, instrument calibrations, and high-resolution geochemical map generation for active fumarolic fields.

Developed at the **Institute of Geophysics, National Autonomous University of Mexico (UNAM)*.

---

## Key Features

- **Calibration Toolbox:** Derive precise sensor coefficients (Gain, Offset, and Cross-Interference) and manage calibration history via an integrated JSON/CSV database.
- **Saturation Correction:** Automatically reconstruct truncated concentration peaks caused by short-term sensor saturation using a linear-extrapolation algorithm.
- **Auto-Lag Correction:** Align non-synchronized sensor time-series through an automated Pearson correlation-based search loop (1–60 s) to determine response delays.
- **Baseline Detrending:** Eliminate environmental or systemic sensor drift using SciPy-based linear detrending while preserving the original signal levels.
- **Cross-Interference Scaling:** Quantitatively remove electrochemical cross-sensitivities (e.g., SO₂ artifacts on H₂S sensors) using real-time mathematical compensation.
- **Advanced Signal Filtering:** Smooth data using multiple options, featuring a optimized **Zero-Phase (Forward-Backward) Single-Pole Low-Pass Filter** to eliminate artificial time lags. Savitzky-Golay and Moving Average filters are also supported.
- **Advanced Linear Regression:** Computes molar gas ratios (R-Squared, MAE) using an ML-optimized framework (80/20 Train/Test split). This approach prevents overfitting and enhances R-Squared coefficient precision.
- **GIS Mapping Utilities:** Synchronize processed gas ratios with GPS data using PCHIP or Linear spatial interpolation to export publication-ready PDF, KML, and GIS-compatible Shapefiles.

---

## Installation & Quick Start

VGR Processor is compiled as a standalone application for Windows, requiring **no prior Python installation** or external library setups.

1. Go to the [Releases](https://github.com/alifhh/VGR-Processor/releases) section of this repository.
2. Download the latest `VGR_Processor.zip` package.
3. Extract the ZIP file into a local folder on your computer.
4. Run the executable file: `VGRP.exe`
   *(Note: The first launch may take a few minutes to initialize the isolated virtual environment).*

---

## Input Data Format

The software accepts data in `.csv` or `.xlsx` (Excel) formats. 
- **Gas Ratios:** Files should contain at least three columns: `Time`, `Gas 1`, and `Gas 2`. Column headers must match standard gas names (e.g., `CO2`, `SO2`, `H2S`, `H2O`) without special characters.
- **GPS Coordinates:** Files should include `Time`, `Longitude`, and `Latitude` columns for proper spatial synchronization.

---

## Documentation

For a comprehensive step-by-step guide on formulas, physical sensor models, and detailed processing workflows, please consult the complete **User Manual** provided in the repository files or attached as Supplementary Material in our manuscript.

---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Citation

If you use VGR Processor in your research or monitoring campaigns, please cite our corresponding paper:
> *Zare, A., Campion, R., Bahmani, A., Paulín Zavala, T. (2026). VGR Processor: A Python-Based Software for Calculating Gas Ratios from Multicomponent Volcanic Gas Analyzers. Computers & Geosciences.*
