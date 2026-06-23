# TNEP-Storage

A Julia package for large-scale **Transmission Network Expansion Planning (TNEP)** with **grid-scale energy storage**, using PTDF-based network formulations and decomposition methods for high-renewable power systems.

## Description 

TNEP-Storage implements optimization models for jointly planning transmission upgrades and energy storage investments under future generation mix and demand scenarios. The framework supports representative-day operational modeling, synthetic Texas test systems, and scalable solution procedures based on continuous relaxations and trust-region Benders decomposition.

The code accompanies the paper:

**High-Resolution PTDF-Based Planning of Storage and Transmission Under High Renewables**

Paper: https://arxiv.org/abs/2510.14696

## Setup

### Installation

1. Clone this repository:

   ```bash
   git clone https://github.com/AI4OPT/TNEP-Storage.git
   cd TNEP-Storage
   ```

2. Install Julia.

3. Create and activate the project environment:

   ```julia
   using Pkg
   Pkg.activate(".")
   Pkg.instantiate()
   ```

This will install all required package dependencies specified in `Project.toml`.

4. To load all source files and dependencies, run:

   ```julia
   include("src_clean/main.jl")
   ```

### Data Setup

Most of the synthetic Texas system data used in this project is obtained from the following dataset:

* https://zenodo.org/records/4538590

Because the dataset is too large to store in this repository, it must be downloaded separately and extracted into the root directory of the project.

After extraction, rename the top-level folder to `tamu`.

The repository structure should look similar to:

```text
TNEP-Storage/
├── src_clean/
├── data/
├── examples/
├── Project.toml
├── README.md
├── ...
└── tamu/
    ├── 2020/
    ├── 2030_ambitious_goals/
    ├── 2030_current_goals/
    ├── base_grid/
```

## Usage

> **Note:** A valid Gurobi license is required to solve the large-scale instances used in this project.

### Configuration Setup

An example configuration file is provided at:

```text
examples/example_simdir/config.toml
```

The most important configuration parameters are shown below:

```toml
dates = ["2016-08-11"]
representative_prob = [1.0]
num_representatives = 1
decarbonization_year = 2022
decarbonization = "data/topology/tamu/decarbonization_1015.csv"
time_series_dir = "tamu/base_grid"
power_system_data = "data/topology/tamu/texas/power_system_data.json"
```

These parameters specify:

* `dates`: Representative days used in the planning model.
* `representative_prob`: Weights associated with each representative day.
* `num_representatives`: Number of representative days.
* `decarbonization_year`: Target year for the generation mix projection.
* `decarbonization`: Path to the projected generation mix data.
* `time_series_dir`: Directory containing load, wind, and solar time-series data.
* `power_system_data`: Network topology and power system data.

Additional configuration options control investment costs, storage parameters, solver settings, trust-region parameters, and other modeling assumptions.

### Warm-Start Procedure

To run the warm-start procedure described in the accompanying paper for a single representative day, execute:

```julia
simdir = "examples/example_simdir"
model, data = run_model(simdir)
```

This command solves the continuous relaxation of the planning model for the representative day specified in the configuration file and returns the corresponding relaxed investment decisions.

The resulting `model` and `data` objects can then be used for further analysis or as inputs to subsequent stages of the solution procedure.

### TODO

Additional functionality and run configurations will be documented in future updates.

## Source Paper

If you use this code or the associated methodology, please cite:

```bibtex
@article{wu2025high,
  title={High-Resolution PTDF-Based Planning of Storage and Transmission Under High Renewables},
  author={Wu, Kevin and Haider, Rabab and Van Hentenryck, Pascal},
  journal={arXiv preprint arXiv:2510.14696},
  year={2025}
}
```

Test


