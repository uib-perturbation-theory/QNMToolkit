This folder includes the following illustrative notebooks:

- **`initial-data.ipynb`**  
  Construction and analysis of localized Gaussian perturbations, including their spectral content and dependence on the parameters $(\sigma, \nu)$.

- **`generic-waveforms.ipynb`**  
  Time-domain numerical evolution of the Regge–Wheeler equation and extraction of waveforms at a fixed observer location, exploring the impact of the initial data on the resulting signal.

- **`theoretical-qnecs.ipynb`**  
  Description of the algorithm used to compute quasi-normal excitation coefficients (QNECs), both in the asymptotic approximation and using the Leaver wavefunction, for the fundamental mode and higher overtones.

- **`numerical-qnecs.ipynb`**  
  Presentation of the numerical procedure used to extract QNECs directly from time-domain waveforms.

- **`waveform-reconstruction.ipynb`**  
  Reconstruction of waveforms using the QNECs derived in `theoretical-qnecs.ipynb`, illustrating the contribution of individual modes to the full signal.

These notebooks form a complete pipeline: initial data → time evolution → QNEC extraction → waveform reconstruction.
