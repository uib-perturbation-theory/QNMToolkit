# Shaping black hole resonances

This repository contains a collection of codes and scripts used for data analysis and analytical derivations associated with the following articles:

> **A. Svyatkovskyy Kholyavka, J. A. León Vega, X. Jiménez Forteza & S. Datta (2026)**  
> *“Shaping black hole resonances I. Black hole ringdown as a spectral filtering process”*  
> [arXiv:gr-qc/2604.XXXXX](https://arxiv.org/abs/gr-qc/2604.XXXXX)

> **A. Svyatkovskyy Kholyavka, J. A. León Vega, X. Jiménez Forteza & S. Datta (2026)**  
> *“Shaping black hole resonances II. Structure and optimization of quasinormal mode excitation”*  
> [arXiv:gr-qc/2604.XXXXX](https://arxiv.org/abs/gr-qc/2604.XXXXX)

## Repository Structure

The repository is organized as follows:

- **`notebooks/`**  
  Contains reproducible Jupyter notebooks illustrating the main ideas and numerical experiments:
  - `initial-data.ipynb`: Construction and spectral analysis of localized Gaussian perturbations, including their dependence on $(\sigma, \nu)$.
  - `generic-waveforms.ipynb`: Time-domain evolution of the Regge–Wheeler equation and extraction of waveforms at a fixed observer location.
  - `theoretical-qnecs.ipynb`: Computation of quasi-normal excitation coefficients (QNECs) using both asymptotic approximations and the Leaver wavefunction, for the fundamental mode and higher overtones.
  - `numerical-qnecs.ipynb`: Extraction of QNECs directly from numerical waveforms.
  - `waveform-reconstruction.ipynb`: Reconstruction of waveforms using the computed QNECs, illustrating the contribution of individual modes.

- **`rw_solver/`**  
  Python implementation of a time-domain solver for the Regge–Wheeler equation.  This module provides the numerical infrastructure used to evolve perturbations, including:
  - Generation of initial data (Gaussian and oscillatory profiles)
  - Finite-difference time evolution
  - Extraction of waveforms at prescribed observer locations  

## Citation

If you use this repository in your research, please cite the corresponding papers:

```bibtex
@article{Svyatkovskyy2026a,
  author        = {Svyatkovskyy Kholyavka, A. and León Vega, J. A. and Jiménez Forteza, X. and Datta, S.},
  title         = {Shaping black hole resonances I. Black hole ringdown as a spectral filtering process},
  year          = {2026},
  eprint        = {2604.XXXXX},
  archivePrefix = {arXiv},
  primaryClass  = {gr-qc},
  month         = {4}
}

@article{Svyatkovskyy2026a,
  author        = {Svyatkovskyy Kholyavka, A. and León Vega, J. A. and Jiménez Forteza, X. and Datta, S.},
  title         = {Shaping black hole resonances II. Structure and optimization of quasinormal mode excitation},
  year          = {2026},
  eprint        = {2604.XXXXX},
  archivePrefix = {arXiv},
  primaryClass  = {gr-qc},
  month         = {4}
}
