# Temperature Diffusion Simulation

This program simulates the diffusion of temperature in a two-dimensional space using the lattice Boltzmann method. The simulation is based on the following parameters:

- `m`: distance where the temperature is diffused
- `time_step`: time step in the simulation
- `space_step`: space step in the simulation
- `viscosite`: the viscosity of the fluid
- `prandlt`: the Prandtl number
- `coeffision`: coefficient of diffusion
- `relaxation_time`: relaxation time of the fluid
- `relaxation_time2`: relaxation time of the temperature
- `cs`: sound velocity
- `constantev0`: constant velocity of the fluid
- `U0`: initial temperature
- `iterations`: number of iterations in the simulation

The simulation consists of two main steps: collision and streaming. In the collision step, the distribution functions of the fluid are updated using the equilibrium distribution functions, while in the streaming step, the distribution functions are shifted to neighboring cells.

Finally, the program plots the two-dimensional space with the diffusion of temperature using the `matplotlib` library.

## Dependencies

This program requires the following Python libraries:

- `numpy`
- `matplotlib`

## How to Use

1. Install the required libraries by running `pip install numpy matplotlib`
2. Copy and paste the code into a Python environment or save it as a `.py` file.
3. Run the program. The diffusion of temperature will be displayed in a plot.