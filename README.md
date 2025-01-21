# Viability_IK
## Introduction
``Viability_IK`` is a quadratic programming (QP) based inverse kinematics (IK) method to simultaneously handle physical joint limits (_including joint range, velocity and acceleration limits_) and avoidance of whole-body collision (_including self-collision and collisions with static obstacles_).

These constraints may conflict with each other and result in infeasible solutions of IK, which can subsequently produce unpredictable motions.

``Viability_IK`` incorporates an additional linear constraint to ensure that the robot state is **viable**, guaranteeing the existence of solutions that do not violate the given constraints.

It can effectively handle most robots with DOFs below 10 and can also accommodate some robots with higher DOFs under simpler constraints.

## Installing within Conda
``Viability_IK`` relies on several optimization solvers, each with complex dependencies, making it difficult to build from source.
Therefore, we strongly recommend using the provided Conda environment for installation, as detailed in the tutorial below.
This approach allows for quick deployment across Windows, macOS, and Linux.
For users who prefer to install everything from source, we have listed the main dependencies at the end.
### Installation of Conda
To get started with conda (or mamba) as package managers, you need to have a base conda installation. Please _do not_ use the Anaconda installer, but rather start with [Miniforge](https://github.com/conda-forge/miniforge) that is much more "minimal" installer. This installer will create a "base" environment that contains the package managers conda and mamba. 

To create a new Conda environment with a specified name (e.g., ``via_ik``) and activate it, run the following commands:
```sh
conda create -n via_ik python=3.11 -y
conda activate via_ik
```
In Conda installations managed by Miniforge, ``conda-forge`` is set as the default (and only) channel. It provides all the necessary dependencies across major platforms:
```sh
conda install scip osqp casadi compilers cmake make git doxygen pkgconfig boost pybind11 matplotlib -c conda-forge
```
### Build and installation in conda
The Build and installation of ``Viability_IK`` in Conda can be performed using the following standard CMake procedure:
```sh
git clone --recursive https://github.com/zyc1155/Viability_IK.git
cd Viability_IK
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release
cmake --build build --target install --config Release
```
## Examples
The ``examples`` folder, located in the root directory, contains two example files: ``example_cxx.cpp`` for C++ and ``example_py.py`` for Python.
### Example for CXX
An executable ``viability_ik_example_cxx``  will be installed into ``${CONDA_PREFIX}/bin`` on Unix-like systems or ``${CONDA_PREFIX}/Library/bin`` on Windows. This executable is generated from   ``example_cxx.cpp``. 
To test the library installation and to verify the IK method, run it from a terminal with the corresponding Conda environment activated.
### Example for PYTHON
A python wrapper of ``Viability_IK`` will be installed into conda env. 
You can launch ``example_py.py`` to have fun with a simulation of two DOF robot.
## Common build
### CMake options
* ``BUILD_IN_CONDA`` Build program in conda enviroment (ON/OFF, default: ON)
* ``BUILD_PYTHON``   Build python wrapper (ON/OFF, default: ON)
* ``BUILD_CXX_EXAMPLE `` Build cxx example (ON/OFF, default: ON)
### Dependencies

