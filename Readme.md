# General Eccentric Binary Pulsar Analysis

This project reproduces and extends results from the paper : 
**Bagchi, Lorimer & Wolfe (2018) — "On the detectability of eccentric binary pulsars" ([arXiv:1302.4914](https://arxiv.org/pdf/1302.4914))**. 

This simulation tool is designed to model the orbital dynamics of a binary pulsar system. It calculates key observable quantities, such as **line-of-sight velocity** and **acceleration**, and computes a signal detectability metric known as `gamma1`, which accounts for the rapidly changing Doppler shift of the pulsar's signal.

This project is structured as a multi-project Visual Studio solution, separating the core scientific library from the specific command-line applications.

##  Methods
- **Kepler’s equation solver**: Newton–Raphson iteration with robust initial guesses.
- **Velocity (Eq. 24)** and **Acceleration (Eq. 26)** from Bagchi et al. (2018)
- **gamma1 function** : Equation 17 in the paper.
- **Parameters configurable**: Pulsar mass, companion mass, orbital period, eccentricity, inclination, longitude of periastron (varpi).
- Fully implemented in C++ 

# Prerequisites
To build and run this project, you will need: 
- Visual Studio 2022 with the "Desktop development with C++" workload installed.
- Git for cloning the repository.


# How to Build the Code
The project is built using the provided Visual Studio solution file, which manages all dependencies and configurations.

1. Clone the Repository
Open a terminal and clone the project to your local machine: 

2. Open the Solution:
Navigate into the project folder and open the PulsarSolution.sln file. Visual Studio 2022 will load the solution with all its projects (libgamma, gamma1_unweighted, etc.).

3. Build the Solution:
From the top menu in Visual Studio, select Build -> Build Solution. You can also use the keyboard shortcut `Ctrl+Shift+B`.

This will compile the libgamma static library first, and then build all the individual applications. The final executable files (.exe) will be placed in the bin/ directory.

# How to Run the Code
Each application is designed to be run from the command line or directory from Visual Studio. They read parameters from a corresponding file in the input/ directory.

## Running from Visual Studio (highly Recommended) 
This is the easiest way to test and run the programs.
1. **CHoose a program to Run**: In the Solution Explorer on the right, right-click on the project you want to run (e.g. gamma1_unweighted)
2. **Set as Startup Project**: From the context menu, select "Set as Startup Project". The projects name will become bold.
3. **Modify Inputs (Optional)**: Before running, you can edit the corresponding file in the input/ folder to change the physical parameters of the simulation.
4. **Run**: Press `Ctrl+F5` (or select `Debug -> Start Without Debugging`) to run the program. A console window will appear showing the program's progress and output.
---
