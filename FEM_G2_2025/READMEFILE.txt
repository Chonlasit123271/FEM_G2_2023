Here is a README file for the provided MATLAB code:

# README

## Overview
This MATLAB code performs a finite element analysis of a 3D structural model of a "Cantilever column with displacement control at the top" using TET4 elements.
It calculates the global mass matrix, stiffness matrix, applies support conditions, computes the displacement and stress/strain distributions, and generates force-displacement graphs and visualizations of the deformed structure.

## Input
The code requires the following input:
- `height`: Height of the structure (m)
- `breadth`: Breadth of the structure (m)
- `width`: Width of the structure (m)
- `nu`: Poisson's ratio
- `E`: Young's modulus of the material (N/m^2)
- `rho`: Density of concrete in kg/m^3
- `g`: Acceleration due to gravity, m/s^2
- `SigmaMax_C` = Maximum compressive stress in N/m^2
- `SigmaMax_T` = Maximum tensile stress in N/m^2  

# --------------------------------------------------------------- #

## Script files:
   1. MeshG2.m (or Mesh file): Mesh file of the model which obtained from GMSH
   2. GetStressHDCompression: Material Model or Stress-Strain curve
   3. G2_Main.m (or Main script): Main script to analyze the non-linear dynamics problem including post-processing. 
   4. MeshG2_2: Applies contour stress plots in the post-processing part.

Let's consider detial of each script files:

1. Mesh file 
- Input mesh file as a matlab file.
- Generate node coordinates, element connectivity, number of elements, number of nodes and number of degree of freedom

2.GetStressHDCompression
- Material Model or Stress-Strain curve

3. Main Script 
3.1 Initialize the element and node data structures.
3.2 Calculate the shape functions, element stiffness matrices, mass matrix, body forces using Gaussian integration.
3.3 Assemble the global stiffness matrix by accumulating the element stiffness contributions.
3.4 Apply support conditions
3.5 Apply the load as an incremental load vector
3.6 To solve this dynamic equation apply newmark beta method 
3.7 Apply gradual lateral displacement at top corner node and update displacement, velocity, and acceleration. 
3.8 Use incremental method with equilibrium correction to account for nonlinearity of the material 
3.9 Calculate the nodal strains and convert it to principle strains
3.10 Find principle stresses using relevent tensile or compressive curve
3.11 Check the failure criteria. If it is failed in tension, reduced the elastic modulus and stiffness of the material
3.12 The code will run until the structure fails due to large deformation
3.13 Save the results to the output files.
3.14 Display the post-processing 
	- Deflected Displacement
	- Display failed elements (animation) and
	- Principal Stresses

4. MeshG2_2: Apply the plot contour stresses in post-processing part.
- Extends the script for plotting contour principal stresses in the post-processing part.
- Can be included in the main script if the user develops the necessary commands.
- To plot contour principal stresses:
	- Element coordinates, connectivity, etc., are used.
	- Results of actual principal stresses (sigma1, sigma2, sigma3) are required.


## Output
The code generates the following output:
- Console:
  - Failure summary (element, increment, applied displacement)
- Figures:
  - 3D plots showing:
    - Undeformed mesh
    - Deformed shape
    - Failed elements (stepwise and cumulative)
- Video:
  - `FailedElementsAnimation.mp4` showing failure propagation over time


## Dependencies
The code uses the following MATLAB functions:
- `lgwt`: Computes the Gauss-Legendre quadrature points and weights.
- `GetHDStressStiffness`: Calculates the stress and stiffness relationship 

## Notes
- The code assumes that the input mesh file follows the specified format.
- The support and load conditions are defined in the code and can be modified as needed.
- The code saves the final results in the MATLAB workspace and generated output files.

Version: G2_15MAY2025

