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

# Running on Windows Machine
## Part A : How to Build the Code
The project is built using the provided Visual Studio solution file, which manages all dependencies and configurations.

1. Clone the Repository
Open a terminal and clone the project to your local machine: 

2. Open the Solution:
Navigate into the project folder and open the PulsarSolution.sln file. Visual Studio 2022 will load the solution with all its projects (libgamma, gamma1_unweighted, etc.).

3. Build the Solution:
From the top menu in Visual Studio, select Build -> Build Solution. You can also use the keyboard shortcut `Ctrl+Shift+B`.

This will compile the libgamma static library first, and then build all the individual applications. The final executable files (.exe) will be placed in the bin/ directory.

## Part B : How to Run the Code
Each application is designed to be run from the command line or directory from Visual Studio. They read parameters from a corresponding file in the input/ directory.

## Running from Visual Studio (highly Recommended) 
This is the easiest way to test and run the programs.
1. **CHoose a program to Run**: In the Solution Explorer on the right, right-click on the project you want to run (e.g. gamma1_unweighted)
2. **Set as Startup Project**: From the context menu, select "Set as Startup Project". The projects name will become bold.
3. **Modify Inputs (Optional)**: Before running, you can edit the corresponding file in the input/ folder to change the physical parameters of the simulation.
4. **Run**: Press `Ctrl+F5` (or select `Debug -> Start Without Debugging`) to run the program. A console window will appear showing the program's progress and output.


# Running on Linux Machine
### Requirements

* A C++17 compiler (e.g. `g++` or `clang++`)
* OpenMP support (usually included with `g++` on Linux)

Check if you have `g++` installed:

```bash
g++ --version
```


## Step 1 : Building

Clone the repository and navigate into it:

```bash
git clone https://github.com/Messiers87/binary-pulsar-analysis.git
cd binary-pulsar-analysis
```

Compile one of the programs, e.g. `v2.cpp` (which is gamma1_weighted code):

You can name the ouput file anything you want ! Run this following code.
```bash
g++ -std=c++17 -fopenmp apps/gamma1_weighted.cpp src/*.cpp -o filename
```


This creates an executable named `filename` in the root directory that is `/PulsarSolution` .

## Step 2 : Running

Run the executable:

```bash
./filename
```

* Input is read from `input/gamma1_weighted.input`
* Results are written to `output/gamma1_weighted.output`


## Output

The program prints progress to the terminal and writes results in tabular form to the corresponding `output/` file.

---
