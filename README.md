# Gmsh2Matlab

## read_gmsh.m

It reads a 2D mesh in .msh file created by Gmsh into the Matlab as a structure **T**. 
Here, **T** has all information about the mesh, and is stated as in chapter 6 
of the book **Understanding and implementing the Finite element method** by 
Mark S. Gockenbach.

Here, I assume the Dirichlet edge has tag 1 and the Neumann edge has tag 2.
