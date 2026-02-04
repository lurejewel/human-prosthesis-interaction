# Predictive Forward Simulation of Human Walking
*A Reflex-Based Musculoskeletal Model Optimized via CMA-ES*
## Overview
This repository provides a predictive, experimental-data-free forward dynamic simulation framework for human walking, based on a muscle-reflex controller and evolutionary optimization.

The framework generates stable, human-like level-ground walking by tuning the parameters of a reflex-based controller using the Covariance Matrix Adaptation Evolution Strategy (CMA-ES). No motion capture data or predefined joint trajectories are required.

The current release (v1.0) demonstrates steady-state walking at 1.0 m/s using a lower-limb musculoskeletal model.


## Key Components

The predictive simulation framework consists of four thigh coupled modules:
- **Musculoskeletal Model**
  - 16 degrees of freedom
  - 14 Hill-type muscles (lower limbs only)
  - Forward dynamic simulation using muscle-driven actuation
- **Gait Phase Detection**
  - Real-time detection of gait phases for each leg
  - Phase information updated continuously from the system state
- **Muscle Reflex Controller**
  - CNS-inspired reflex mechanisms
  - Muscle excitations generated from proprioceptive feedback and gait phase
  - Reflex parameters shared between left and right homologous muscles
- **Optimization via CMA-ES**
  - Reflex parameters optimized to achieve stable and efficient walking
  - Parallelize evaluation of candidate solutions
  - Fitness based on gait stability and performance metrics


## Code Structure
```
.
├── demo_predictiveForwardSimulation_humanModel.m   % Main entry script
├── assets/                                        % Data and auxiliary files
├── model/                                         % Musculoskeletal model files
├── functions/                                     % Simulation, control, and optimization functions
└── README.md
```

The main script `demo_predictiveForwardSimulation_humanModel.m` performs:
1. Environment setup and parallel pool initialization
2. Simulation and model configuration
3. CMA-ES optimizer initialization
4. Iterative loop of :
   - Parameter sampling
   - Forward dynamic simulation
   - Fitness evaluation
   - CMA-ES parameter update

## Requirements
### Software
- **MATLAB R2024b**
- Parallel Computing Toolbox
- **OpenSim 4.0 or later**
  - OpenSim must be installed and its **MATLAB API properly configured**
> **OpenSim-MATLAB Interface
> This project relies on the OpenSim MATLAB API for musculoskeletal modeling and forward dynamic simulation. Instructios for configuring the OpenSim-MATLAB interface can be found at: https://opensimconfluence.atlassian.net/wiki/spaces/OpenSim/pages/53089380/Scripting+with+Matlab
> After successful configuration, OpenSim classes (e.g., `org.opensim.modeling.*`) should be accessible directly from MATLAB.

### How to Run
1. Clone this repository and open it in MATLAB.
2. Start the simulation by running:
```
demo_predictiveForwardSimulation_humanModel
```
3. The optimizer will iteratively search for reflex parameter sthat produce stable walking at the target speed.

Optimization results are generated automatically based on the configuration.

## Simulation and Optimization Settings
Key parameters defined in the main script include:
- **Target walking speed**: 1.0 m/s
- **Simulation duration**: 10 s
- **Integration time step**: 0.005 s
- **Parallel workers**: configurable (default: 6)

The optimization terminates when:
- Fitness improvement stagnates over a predefined number of generations (100), or
- The maximum number of generations (1000) is reached.

Due to the stochastic nature of CMA-ES, optimization results may vary slightly across runs.

## Scope and Limitations
- Focused on level-ground steady-state walking
- Reflex parameters are identical for left and right muscles
- Controller structure and parameter bounds are tuned specifically for the current musculoskeletal model (H0714)
- Upper-body dynamics are not included

## Planned Extensions
- **Independent reflex parameters for bilateral muscles**
(next minor release)
- **Extended musculoskeletal models with additional muscle actuators**
(subsequent minor releases)
- **Human-prosthesis coupled models**
(next major release)
- **Additional motor tasks beyond walking**
*(e.g., slope walking, running; release timeline to be determined)*

## Citation
If you use this code in your research, please cite the corresponding publincation:

Jin W, Liu J, Zhang Q, et al. Forward dynamics simulation of a simplified neuromuscular-skeletal-exoskeletal model based on the CMA-ES optimization algorithm: framework and case studies[J]. Multibody System Dynamics, 2024, 62(4): 525-558.

## Contact
Wei JIN

Peking University

wjin24@stu.pku.edu.cn