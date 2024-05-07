# Finite Element Analysis for 3D Tetrahedral Elements


This MATLAB script performs a Finite Element Analysis (FEA) on 3D tetrahedral elements. 
It computes the global stiffness matrix, applies boundary conditions, and solves for node displacements under linear incremental displacement loading (Can be linear cyclic or monotonic load).
Additionally, it includes visualization of both undeformed and deformed structures in both matlab and paraview.

## Requirements

- MATLAB R202x (where x is your version)
- Gmsh for mesh generation (".m" file should be compatible)
- Excel file ".xlsx" file for displacement loading data

## How to Run

1. Prepare your mesh in Gmsh and export it as a " .m" file
2. Manually create the displacement load increment file in " .xlsx" file 
3. Ensure the " .m" and ".xlsx"  are in the same directory as the script.
4. Run the script in MATLAB. Make sure all dependent files are in MATLAB's current directory or in a path recognized by MATLAB.

## Code Structure

- **Initialization**: Clears variables and figures, and initializes the timer.
- **Input Parameters**: Material properties such as Young's Modulus and Poisson's Ratio are defined.
- **Mesh Import**: The script imports node coordinates and element connectivity from a Gmsh-generated " .m" file.
- **Supports and Loads**: Defines boundary conditions and applied loads.
- **Stiffness Matrix Calculation**: Constructs the local and global stiffness matrices for the structure.
- **Solver**: Applies loads and computes the displacements using the global stiffness matrix.
- **Visualization**: Plots the original and deformed geometries to visualize the displacement.

## Outputs

- Plots of undeformed (blue) and deformed (red) structures.
- Displacement vectors for each node.
- Stress calculation at specified points within each element.
- A ".vtk" file is exported to be open the displacements and stresses visualization in Paraview. 

## More Detail can be seen in the script

##The script was developed by First year master's student from Asian Institute of Technology(AIT) for Finite Element Method for Engineering Project