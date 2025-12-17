# Wave-Propelled USV Model: AutoNaut Simulation

**Project for NTNU Course: TTK4550 Engineering Cybernetics, Specialization Project**

This repository contains a MATLAB implementation of a **wave-propelled uncrewed surface vehicle (USV) model**, specifically for the **AutoNaut**. The model is based on the work of [Tufte (2025)](https://torarnj.folk.ntnu.no/CAMS_2025_Wave_Propelled_USVs_FINAL.pdf), developed as part of a specialization project at NTNU.

---

## Overview
The AutoNaut is a wave-propelled vessel designed for long-duration ocean monitoring. This project implements Tufte’s model, combining **traditional ship manoeuvring theory** with **foil theory** to simulate the vessel’s dynamics. The model is designed for ease of use in further research, control system development, and mission planning.

---

## Current Status (17th December 2025)
- The **foil dynamics** in the current implementation **do not work properly** and require further development and validation.

---

## Features
- **Hybrid Model**: Combines ship manoeuvring and foil theory for accurate simulation of wave-propelled USVs.
- **MATLAB Implementation**: Ready-to-use scripts for simulating vessel dynamics.
- **Modular Design**: Easy to extend for control systems, validation, or integration with other tools.

---

## Dependencies
- **MATLAB** (R2020b or later recommended)
- **[MSS Toolbox](https://github.com/cybergalactic/MSS)**: Required for marine systems simulations.

Install the MSS Toolbox by cloning the repository and adding it to your MATLAB path:
```matlab
addpath('path/to/MSS');
```

---

## Usage
1. Clone this repository:
```bash
git clone https://github.com/tfossdal/Autonaut.git
```
2. Open MATLAB and navigate to the project directory.
3. Run the main simulation script:
```matlab
SIMautonaut.m
```
4. Change the parametres of the vessel at the begging of the autonaut.m file
---

## Acknowledgments
* **Tufte, A.** for doing all the modelling, and helping me the maths and physics thoughout the project.
* **Gryte, K.** for being my advisor and helping me complete the project.
