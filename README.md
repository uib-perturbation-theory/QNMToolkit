# Shaping black hole resonances

This repository contains a collection of codes and scripts used for data analysis and analytical derivations associated with the research article:

> **A. Svyatkovskyy Kholyavka, J. A. León Vega, X. Jiménez Forteza & S. Datta (2026)**  
> *“Shaping black hole resonances I. Black hole ringdown as a spectral filtering process”*  
> [arXiv:gr-qc/2605.XXXXX](https://arxiv.org/abs/gr-qc/2605.XXXXX)

## Repository Structure

The repository is organized as follows:

- **`notebooks/`**  
  Contains reproducible Jupyter notebooks illustrating the main ideas and numerical experiments:
  - **`initial-data.ipynb`**  
  Construction and analysis of localized Gaussian perturbations, including their spectral content and dependence on the parameters $(\sigma, \nu)$.
- **`generic-waveforms.ipynb`**  
  Time-domain numerical evolution of the Regge–Wheeler equation.
- **`theoretical-qnecs.ipynb`**  
  QNECs both in the asymptotic approximation and using the Leaver wavefunction, for the fundamental mode and higher overtones.
- **`numerical-qnecs.ipynb`**  
  Numerical procedure used to extract QNECs directly from time-domain waveforms.

- **`rw_solver/`**  
  Python implementation of a time-domain solver for the Regge–Wheeler equation.  This module provides the numerical infrastructure used to evolve perturbations, including:
  - Generation of initial data (Gaussian and oscillatory profiles)
  - Finite-difference time evolution
  - Extraction of waveforms at prescribed observer locations  

## Citation

If you use this repository in your research, please cite the corresponding paper:

```bibtex
@article{Svyatkovskyy2026,
  author        = {Svyatkovskyy Kholyavka, A. and León Vega, J. A. and Jiménez Forteza, X. and Datta, S.},
  title         = {Shaping black hole resonances I. Black hole ringdown as a spectral filtering process},
  year          = {2026},
  eprint        = {2605.XXXXX},
  archivePrefix = {arXiv},
  primaryClass  = {gr-qc},
  month         = {5}
}
