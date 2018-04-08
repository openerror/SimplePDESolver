## Simple PDE solve for reaction-diffusion systems
The forward-Euler method is implemented for solving reaction-diffusion PDEs; in principle, an arbitrary number of species can be accommodated. The initial commit has been tested on simple diffusion, and the non-linear system described in [1].

## TODO
1. Optimize of the Laplacian operator, by avoiding the use of if/else statements

## References
1. Haselwandter, Christoph A., et al. “Self-Assembly and Plasticity of Synaptic Domains through a Reaction-Diffusion Mechanism.” Physical Review E, vol. 92, no. 3, Sept. 2015, p. 032705. APS, doi:10.1103/PhysRevE.92.032705.
