# QNMToolkit

QNMToolkit is a numerical algorithm for studying the quasi-normal excitation coefficients (QNECs) of black hole ringdown signals. The code included in this repository was developed in the context of the research article:

> **A. Svyatkovskyy Kholyavka, J. A. León Vega, X. Jiménez Forteza & S. Datta (2026)**  
> *“Shaping black hole resonances I. Black hole ringdown as a spectral filtering process”*  
> [arXiv:gr-qc/2605.XXXXX](https://arxiv.org/abs/gr-qc/2605.24704)

## Repository structure

The repository is organized as follows:

- **`numerical-qnecs.ipynb`**  
  Jupyter notebook containing the QNMToolkit algorithm used to extract QNECs from black hole ringdown waveforms.

- **`rw_solver/`**  
  Python implementation of a time-domain solver for the Regge–Wheeler equation, used to generate the ringdown waveforms analyzed in the notebook.

## Requirements

The notebook requires standard Python scientific libraries:

```bash
pip install numpy scipy matplotlib jupyter
```

## Questions and contributions

For questions, suggestions, or contributions, feel free to contact the authors or open an issue in the repository.

## Citation

If you use this repository in your research, please cite the corresponding paper:

```bibtex
@article{Svyatkovskyy2026,
  author        = {Svyatkovskyy Kholyavka, A. and León Vega, J. A. and Jiménez Forteza, X. and Datta, S.},
  title         = {Shaping black hole resonances I. Black hole ringdown as a spectral filtering process},
  year          = {2026},
  eprint        = {2605.24704},
  archivePrefix = {arXiv},
  primaryClass  = {gr-qc},
  month         = {5}
}
