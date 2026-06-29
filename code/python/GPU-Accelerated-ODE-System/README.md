# GPU Accelerated Sensitivity Analysis

Initial experiments with sensitivity analysis resulted in hours-long computations. For a sample size of $2^{20}$, initial simulations took up to 5 hours to complete. The next step was to parallelize the simulations on CPU, and disable unnecessary second-order calculations, which reduced computations for $2^{20}$ samples down to 27 minutes. Using JAX and Diffrax, a GPU vectorized ODE solver was implemented to solve ODE systems for a batch of parameter values. This has brought computations down to just __.

## Installation

To run the Python code, you will need to install the necessary packages. Python 3.10 or higher is recommended. The code was developed for and has only been tested on Linux systems with NVIDIA graphics cards supporting CUDA 12.

```bash
git clone https://github.com/Winthrop-MATH-Alzheimers-ADAPT-SC26/Alzheimers-Research
cd Alzheimers-Research
python3 -m venv venv
source venv/bin/activate
pip install -r requirements-gpu.txt
```