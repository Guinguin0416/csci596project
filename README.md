# Final Project: Visualizing Simulations (Extending Assignment 7)

## Task 3: Animate parallel MD code `pmd.c`

### Introduction
This repository contains the source code and related files for a parallel molecular dynamics (MD) simulation of Lennard-Jones systems using the Message Passing Interface (MPI) standard. The primary objective of this project is to simulate the dynamic behavior of particles under the influence of Lennard-Jones potential and to analyze the system's evolution over time.

### Major Methods and Implementations
The simulation is implemented in C, leveraging the MPI library for parallel computation. The main components of the project include:
- `pmd.c`: The main source file containing the logic for the MD simulation.
- `pmd.h`: Header file defining constants, global variables, and function prototypes.
- `pmd.in`: Input file specifying initial parameters for the simulation.
- `pmd.sl`: A script for submitting the simulation job to a high-performance computing cluster.

#### Overview
The pmd.c program simulates the dynamics of particles interacting via the Lennard-Jones potential using parallel computing with MPI (Message Passing Interface). It models a system of particles in a three-dimensional space where each particle's movement is influenced by its interactions with other particles.

#### Input
The program takes its initial configuration from the `pmd.in` file. This file contains parameters such as:
- The initial unit cell dimensions (InitUcell).
- Density (Density).
- Initial temperature (InitTemp).
- Time step size (DeltaT).
- Total number of time steps to simulate (StepLimit).
- Frequency of reporting average properties (StepAvg).
These parameters are read in the init_params() function.

#### Initialization
- **Setting up the simulation box:** The program initializes the positions of atoms in a face-centered cubic (FCC) lattice and assigns random velocities based on the specified temperature.
- **Domain decomposition:** The simulation space is divided among MPI processes. Each process handles a subset of the entire system.

#### Main Simulation Loop
- **Time stepping:** The single_step() function performs the core of the simulation, updating the positions and velocities of the atoms at each time step using the Velocity-Verlet integration method.
- **Inter-process communication:** Atoms that move across the boundaries of a process's domain are transferred appropriately to neighboring processes (atom_move()).
- **Force computation:** The compute_accel() function calculates the forces acting on each atom due to its neighbors within a certain cutoff radius (defined by the Lennard-Jones potential).

#### Output
- **Data Writing:** The write_xyz() function is called to write the coordinates of the atoms to an output file (output.xyz). This function is invoked by the MPI process with rank 0, ensuring that only one process writes to the file. The output format is suitable for visualization tools to animate the particle movements over time.
- **Property Evaluation:** The eval_props() function evaluates and prints physical properties like temperature, potential energy, and kinetic energy at specified intervals (StepAvg).

#### Parallel Computing Aspects
- **MPI Initialization:** The program begins with initializing MPI and determining the rank of each process.
- **Data Distribution:** Each MPI process works on a subset of the entire atomistic system, which involves calculations specific to its domain and communication with adjacent domains.
- **Synchronization:** MPI barriers and send/receive operations ensure proper synchronization and data transfer among processes.

#### Closing the Simulation
- After completing the specified number of time steps, the program finalizes the MPI environment and exits. The output.xyz file contains the trajectory data which can be used to generate animations or further analyze the molecular dynamics.


### Key Features
- Efficient parallel computation across multiple processors.
- Implementation of the velocity-Verlet integration scheme for time evolution.
- Periodic boundary conditions and spatial decomposition for handling large numbers of particles.

### Results
The results of the simulation are output in the `output.xyz` file, which records the positions of particles at each timestep. This data can be used to visualize the particle trajectories and analyze the system's behavior over time.

### Animation
An animation of the simulation is provided in the repository (in `.mov` format). This animation illustrates the dynamic evolution of the particle system.

<img src="/Task3_AnimatePMD/pngs/untitled.00001.png" width="300" height="300">

### Discussion
The animation of the `pmd.c` simulation offers critical insights into the molecular dynamics of particles interacting via the Lennard-Jones potential. This visualization is instrumental in understanding several key aspects:
- **Particle Dynamics and Interactions:** The animation vividly illustrates how particles move and interact over time. Watching the particles' trajectories and behaviors provides an intuitive understanding of the physical principles at play, such as attraction and repulsion governed by the Lennard-Jones potential.
- **Temporal Evolution:** The animation allows for the observation of the system's evolution over time. This is crucial for identifying transient phenomena and understanding how the system reaches equilibrium, or how it responds to changes in initial conditions.
- **Impact of Simulation Parameters:** The animation makes it easier to comprehend the impact of various simulation parameters, such as particle density, initial temperature, and time step size. These factors critically influence the behavior and stability of the simulated system.
- **Enhancing Comprehension:** As demonstrated in the animation, visual representations are immensely powerful in enhancing our comprehension of complex systems. They allow us to observe and analyze patterns and behaviors that might be non-trivial to discern from raw data or numerical outputs.

In conclusion, the animation of the `pmd.c` simulation is not just a visual aid but a powerful tool for scientific inquiry and education, providing a window into the microscopic world of molecular dynamics.

### How to Run
1. Compile the program using an MPI-compatible compiler:
   ```sh
   mpicc -O -o pmd pmd.c -lm
2. Run the program with slurm script: `pmd.sl`
    ```sh
    sbatch pmd.sl

### Requirements

- An MPI implementation (e.g., OpenMPI, MPICH).
- A C compiler (e.g., GCC).
