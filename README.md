# Particle In Cell - Electrostatic model

Basic (1D) simulator modeling evolution of systems of charged particles (electrostatic) over time. For example, plasma simulations. Uses the particle in cell (PIC) approach to model system state and compute resulting dynamics.

## Model assumptions

Uses uniform Rho, Dirichlet boundary conditions (0V in this case) for electric field model.


## Sources

This is a rust implementation of the earliest model from the book "Plasma Simulations by Example" by Brieda.


## Execution

To run:
```bash
cargo run
