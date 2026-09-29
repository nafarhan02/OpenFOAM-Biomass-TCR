# Biomass TCR Solver for OpenFOAM v2412

Custom OpenFOAM solver developed for biomass combustion modelling.

The solver is based on `coalChemistryFoam` and has been modified to use a Two-Competing-Rate (TCR) devolatilisation model.

## Main modifications

- Custom biomass particle model
- Custom biomass cloud
- Custom `biomassChemistryFoam` solver
- Two-Competing-Rate devolatilisation
- Particle reaction progress tracking using `zeta`
- Temperature-dependent volatile release
- Transition from devolatilisation to char oxidation

## Why TCR?

Biomass particles experience different temperatures during combustion.

The TCR model uses two competing reactions instead of one reaction rate, allowing the devolatilisation behaviour to respond better to changes in particle temperature and thermal history.

## Current status

The solver has been tested using simplified cases.

The customised solver was first compared with the standard OpenFOAM solver under the same conditions.

The TCR model was then tested under low and high temperature conditions to confirm that the model responds as expected.

This is currently solver verification and functional testing only.

The solver has not yet been validated against experimental boiler data or applied to the final boiler geometry.

## Future development

The next stage is to improve the particle model for larger palm biomass fibres and particles, which behave differently from the smaller near-spherical coal particles normally used in OpenFOAM.

## OpenFOAM version

OpenFOAM v2412

## Compilation

Compile the custom library:

    wmake libso

Then compile the solver:

    cd solver/biomassChemistryFoam
    wmake

