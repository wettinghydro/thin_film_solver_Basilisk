# An efficient and highly scalable solver for modelling thin film flows in heterogeneous environments

This repository contains the source code and data generation scripts used to reproduce the figures and results presented in our paper. All code is written in C and uses the Basilisk solver.

## Repository Structure

* **`Figure_XX/`**: Contains the specific code used to create the figure.
* **`case/`**: Contains all the source code for the simulations.
* **`data/`**: This directory is where the simulation data is stored. Due to large file sizes, the raw data is not included in the repository but can be regenerated using the provided scripts.

## Case file explanation
The main files are:

- ```compile.sh```:                bash script to compile the code.
- ```launch.sh```:                 bash script to compile the code with the appropriate flags.
- ```thinFilm.c```:                main source code file.
- ```utils```:                     directory containing auxiliary header files.
- ```parameters.dat```:            input file generated using [utils/inputDataBuilder.c](utils/inputDataBuilder.c).

---

## Getting Started

### Prerequisites
To run the simulations, you need to have the [Basilisk solver](http://basilisk.fr/) installed.

### Compilation
The simulation code can be compiled using the provided `compile.sh` script.

`./compile.sh thinFilm.c + [Compilation Flags]`

For example, to compile with OpenMP (activated by default):

`./compile.sh thinFilm.c`

To compile with MPI:

`./compile.sh thinFilm.c -D_MPI=1`

### Compilation Flags
Additional flags can be used to customize the compilation:

* `-DINCLINED`: Includes gravity effects.
* `-DPRINT`: Outputs flow field data in a Py-compatible format.
* `-DBDF2`: Uses a 2nd-order backward differencing time discretization.

---

## Reproducing Figures

Each subdirectory corresponds to a specific figure in the paper. To reproduce the data and plot for a figure:

1.  Navigate to the specific figure directory (e.g., `Figure_01/`).
2.  Run the `launch.sh` script located in the `case/` directory to generate the required data.
3.  Compile the LaTeX `.tex` file.

For all cases, the post-processed data are included to be used directly in the `.tex` file.

---

## Citation
If you use this code or the concepts from the paper in your own work, please cite the following publication:

```bibtex
﻿@article{Vrionis_2025_ThinFilmSolver,
author={Vrionis, Panayiotis Yiannis
and Demou, Andreas D.
and Karapetsas, George
and Papageorgiou, Demetrios T.
and Savva, Nikos},
title={An efficient and highly scalable solver for modelling thin film flows in heterogeneous environments},
journal={Theoretical and Computational Fluid Dynamics},
year={2025},
month={Sep},
day={13},
volume={39},
number={5},
pages={37},
abstract={We present an open-source, modular and highly scalable thin-film flow solver that can perform complex simulations of wetting and dewetting phenomena in the presence of substrate or environmental heterogeneities. The implementation is based on the thin-film approximation and circumvents the stress singularity at moving contact lines by assuming the presence of an ultra thin precursor film covering the whole surface. The solver is implemented within the open-source Basilisk library (http://basilisk.fr/), and is validated by comparing with the predictions of reduced-order models derived from matched asymptotics analyses for a number of representative cases with isolated droplets. The capabilities of the solver are demonstrated by simulating multiple interacting droplets sliding on an inclined and chemically heterogeneous substrate.},
issn={1432-2250},
doi={10.1007/s00162-025-00757-x},
url={https://doi.org/10.1007/s00162-025-00757-x}
}


```
