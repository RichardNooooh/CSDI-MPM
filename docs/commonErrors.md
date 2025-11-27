# Common Errors and Problems

## `Vertex found outside of the grid` Error

This error appears due to two primary reasons:
1. The initial mesh is not within the bounds of the background grid at initialization
    - Fix: check the grid/mesh dimensions and sizes, then also check if the mesh is in the first quadrant/octant.
2. The mesh numerically explodes.
    - Fix: check the (approximate) CFL condition, modify the material parameters, and/or lower the timestep.

## "Imaginary Number" Error

This will likely arise inside of the constitutive equation of your choice and the errors that appear due to a negative Jacobian, $J = \det(\mathbf{F})$. 

This may occur due to the mesh numerically exploding in some fashion before the "Vertex found outside of the grid" error appears, or it may appear simply due to what you're modeling.

You may benefit from a finer mesh (either the grid and/or the material).

## Pressure oscillations or volumetric locking

This is a common phenomena in explicit timestep MPM due to high stiffness and/or the null space error. You'll have to relax the material parameters to mitigate it or modify the "points per cell" ratio (see paper).
