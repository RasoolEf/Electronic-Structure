DFT Calculations for NiO Systems
This repository contains a set of scripts and configuration files for performing Density Functional Theory (DFT) calculations on NiO systems. The calculations include band structure, density of states (DOS), d-band center, and orbital-projected density of states (PDOS). The repository is designed to work with the GPAW DFT code.

Files Included
1. dband_calc.py
This Python script is used for calculating the band structure and d-band center of the system. It provides essential information about the electronic properties of NiO surfaces.

2. pdos.py
This script calculates the orbital-projected density of states (PDOS) for the NiO systems. It helps analyze the contribution of different orbitals (e.g., Ni 3d, O 2p) to the total density of states.

3. input_config.json
This JSON file contains the configuration settings for generating the NiO systems. The file is designed to create the following NiO surfaces:

NiO(111)

NiO(100)

NiO(110)

4. gpaw_config.json
This JSON file includes the configuration for setting up DFT calculations using the GPAW code. It contains parameters necessary to run the calculations such as k-point mesh, functional, and other simulation settings.

How to Use
Prepare Your System:

Modify the input_config.json file to specify the desired NiO surface (NiO(111), NiO(100), or NiO(110)).

Configure DFT Calculations:

Edit the gpaw_config.json to set up the desired parameters for your DFT calculations, such as the k-point grid, exchange-correlation functional, and other simulation details.

Run DFT Calculations:

Run the DFT calculations using the GPAW package with the provided configurations.

Post-Processing:

Use dband_calc.py to analyze the band structure and d-band center.

Use pdos.py to compute and visualize the orbital-projected density of states (PDOS).

Requirements
GPAW DFT software

Python 3.x

Dependencies: numpy, matplotlib, ase, gpaw

License
This project is licensed under the MIT License.

