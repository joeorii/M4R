# Jacquet–Langlands Computations

This repository contains Python code used for explicit computations related to quaternionic modular forms and the Jacquet–Langlands correspondence, as discussed in Section 7 of the thesis.

## Files

* `JL_algorithm.py` contains the main algorithm.
* `example_d11.py` runs the example with discriminant 11 and level 1.
* `example_d19.py` runs the example with discriminant 19 and level 1.
* `environment.yml` specifies the dependencies required to run the code.

## Requirements

The code requires:

* Conda
* Python 3
* SageMath
* NumPy
* pandas

The required dependencies are specified in `environment.yml`.

SageMath is not supported as a native Python package on Windows. The code is therefore intended to run on macOS or Linux. Windows users should use the Windows Subsystem for Linux (WSL).

## Installation

### 1. Clone the repository

Open a terminal and run:

```bash
git clone https://github.com/joeorii/M4R.git
cd M4R
```

### 2. Install Conda

Install a Conda distribution such as Miniconda, Anaconda, or Miniforge.

Instructions for installing Conda and managing environments are available in the [Conda documentation](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html).

### 3. Create the Conda environment

From inside the repository directory, run:

```bash
conda env create -f environment.yml
```

This creates the Conda environment named `m4r`.

## Activating the environment

To activate the Conda environment, run:

```bash
conda activate m4r
```

## Usage

To run the example with discriminant 11 and level 1, use:

```bash
python example_d11.py
```

To run the example with discriminant 19 and level 1, use:

```bash
python example_d19.py
```

## Deactivating the environment

When finished, leave the Conda environment by running:

```bash
conda deactivate
```
