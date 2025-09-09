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
@article{Vrionis_2025_ThinFilmSolver,
  author  = {P.-Y. Vrionis and Andreas D. Demou and George Karapetsas and Demetrios T. Papageorgiou and Nikos Savva}
  title   = {An efficient and highly scalable solver for modelling thin film flows in heterogeneous environments},
  journal = {Theoretical and Computational Fluid Dynamics},
  year    = {2025},
  doi     = {10.1007/s00162-025-00757-x},
}
```
