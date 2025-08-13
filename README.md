# Faceted Nanoparticle Self-Assembly Simulation using Analytical van der Waals Potentials

This repository provides a Monte Carlo simulation framework for modeling the self-assembly of faceted nanoparticles (NPs) using an accurate **analytical potential** to capture orientation-dependent van der Waals (vdW) interactions. The framework enables large-scale and fast simulations of NP assembly, supporting multiple particle shapes and interaction models.

---

## Key Features

- Accurate closed-form **analytical vdW potentials** between faceted NPs.
- Support for multiple particle shapes: cubes, rods, disks, and triangular prisms.
- Monte Carlo simulation using both **single-particle** and **virtual cluster moves** for improved sampling.
- Comparison with atomistic and coarse-grained (CG) interaction models.
- Capability to simulate **phase behavior** and assembly morphology under varying interaction strengths.

---

## ⚙️ Compilation

Make sure you have `g++` installed (GCC for C++).

To compile the code, run:

```bash
g++ -o my_program vdwmontecarlo.cpp PolygonClipping.cpp Polygon.cpp potential_cal.cpp
```

## ▶️ Running the Simulation
Once compiled, run the simulation with:

```bash
./my_program
```

The simulation starts based on the parameters specified in the main() function of vdwmontecarlo.cpp.

## 🧠 What You Need to Modify

To customize your simulation, you only need to modify the variables in the `main()` function of `vdwmontecarlo.cpp`. Below is a list of configurable variables and their purposes:

- **`Style`** – Sets the visualization style for trajectory output (e.g., `"Virtual"`, `"Vertex"`, `"AA"`).

- **`model`** – Chooses the energy calculation model:  
  - `"vdW"` for analytical van der Waals potential  
  - `"AACG"` for atomistic or coarse-grained models

- **`NumberOfParticles`** – Specifies how many nanoparticles are included in the simulation.

- **`BoxLength`** – Defines the length of the cubic simulation box.

- **`Temperature`** – Sets the system temperature, used to compute thermal energy (`kBT`).

- **`globalminEnergy`** – Defines the desired energy minimum between two particles, controlling interaction strength.

- **`AASigma`** – Diameter of an atom in the atomistic model (used to define particle size).

- **`AAEpsilon`** – Depth of the Lennard-Jones potential well (controls interaction strength between atoms).

- **`Hamaker`** – Hamaker constant used in the analytical van der Waals model (controls interaction strength).

- **`AtomDensity`** – Atomic number density used in the vdW energy calculation (atoms per Å³).

- **`MaxRotation`** – Maximum rotation angle allowed per Monte Carlo trial move (in radians).

- **`MaxMoveSingle`** – Maximum translation allowed for a **single-particle** move.

- **`MaxMoveCluster`** – Maximum translation allowed for a **cluster** move (multiple particles together).

- **`SingleMoveProbability`** – Probability of choosing a single-particle move instead of a cluster move.

- **`EquilibrationSteps`** – Number of initial Monte Carlo steps for equilibration (no data collected during this phase).

- **`ProductionSteps`** – Number of steps during which simulation data is recorded.

- **`EnergyOutputInterval`** – Interval (in steps) at which energy values are printed to the output.

- **`TrajectoryInterval`** – Interval (in steps) at which trajectory data is saved to file.

- **`RestartFileInterval`** – Interval (in steps) for saving restart files (used to resume simulations later).

- **`volumeMoveInterval`** – Interval (in steps) at which volume expansion moves are attempted (if enabled).

- **`CubeSideAASize`** – Number of atoms per side in a particle’s cubic lattice.

- **`aspectratio`** – Height-to-length ratio of a particle (e.g., >1 = rod-like, 1 = cube).

- **`volumeexpansion`** – Enables/disables volume expansion moves.

## 📁 Output Files

The simulation generates:

- **`.lammpstrj` trajectory files** (for visualization)
- **Energy logs** recorded every `EnergyOutputInterval` steps
- **Optional restart files** saved every `RestartFileInterval` steps

<img width="226" height="218" alt="image" src="https://github.com/user-attachments/assets/3b4850c7-eee2-4b5c-b4c1-2cd7998f52f1" />

## 📝 Dependencies
No external libraries required. All geometry and clipping operations are implemented internally (e.g., using Greiner–Hormann algorithm for polygon clipping).


