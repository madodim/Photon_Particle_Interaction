# Photon_Particle_Interaction
Numerical simulation of a photon beam interacting with particles in a medium, using Monte Carlo methods

This project simulates the transport of 10⁶ photons in different geometries (a 2D circle, a 2D infinite rectangle and a 3D sphere) using a Monte Carlo methods. The photons undergo scattering events and can either escape, remain trapped inside the medium, or (in the case of the 2D circle) be absorbed based on an assigned probability.

# Features
📌 Random Number Check: Before running the simulation, the code verifies the behavior of the "rand" function by computing the mean value, standard deviation, and plotting the distribution of random numbers to ensure uniformity <br />
📌 Circle Simulation: Photons start at the center of a circle and scatter until they escape or reach the maximum number of scatterings, and may be absorbed based on a given probability. <br />
📌 Sphere Simulation: Photons start at the center of a sphere and scatter until they escape or reach the maximum number of scatterings. <br />
📌 Rectangle Simulation: Photons enter from the bottom surface (z = 0), scatter inside the medium, and can either escape from the upper surface (z = D) or be reflected towards the bottom surface. "D" is the thickness of the rectangle. <br />
📌 Visualization: Plots show the trajectories of escaped and trapped photons. <br />
📌 Statistical Analysis: The probability of photon escape, trapping, absorption, and scattering is calculated and displayed. <br />

# Requirements
- MATLAB (Recommended version: R2021a or later)
- Basic MATLAB toolboxes

# Results
The simulation outputs:
- The percentage of photons that escape, are absorbed, get trapped, or scatter before escaping.
- 2D and 3D plots visualizing photon trajectories.

# Changes
You can modify the simulation parameters to explore different scenarios:
- different number of photons (increase for better statistics, decrease for faster execution)
- different maximum number of scatterings (works like a time step)
- different optical depth (adjusts mean free path)
- different absorption probability
- different geometry size: Modify the radius (R) of the circle/sphere or the thickness (D) of the rectangle to test different configurations.
