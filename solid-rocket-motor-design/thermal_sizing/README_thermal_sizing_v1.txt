================================================================================
PROJECT: COOLING SYSTEM SIZING
================================================================================

This package contains the MATLAB code for the sizing and optimisation of the 
regenerative cooling jacket and the Thermal Barrier Coating (TBC)

--------------------------------------------------------------------------------
SOFTWARE COMPATIBILITY & PREREQUISITES
--------------------------------------------------------------------------------
Ensure you are using the same version or a compatible one for optimal performance 
and to avoid potential syntax errors.

The code `thermal_sizing_v1.m` has been developed and tested using:
- MATLAB R2025b

--------------------------------------------------------------------------------
HOW TO RUN THE CODE `thermal_sizing_v1.m`
--------------------------------------------------------------------------------
1. DOWNLOAD THE FILE: Ensure you have downloaded the file (.m) to your computer.

2. OPEN MATLAB: Launch MATLAB R2025b.

3. Open the file `thermal_sizing_v1.m`.

4. RUN THE MAIN SCRIPT:
   You can run it by clicking the "Run" button (green arrow) in the "Editor" 
   tab. Select 'add to path' if asked.

5. RUN THE CODE FROM THE TERMINAL:
   - set your current directory to the one containing the code using the command 'cd' (e.g. cd '/Users/namesurname/.../Propulsion/project/');
   - type `thermal_sizing_v1` in the Command Window 
   - press Enter.

   NOTE: Do not attempt to run the functions individually. Only the main script is designed for direct 
   execution;
   The computation may lasts up to one minute.

--------------------------------------------------------------------------------
RESULTS
--------------------------------------------------------------------------------
After running the file `thermal_sizing_v1.m` you will be able to read the results computed. The results consists in two figures,`TBC justification` and `Design Search Space`, and a report printed in the command window. 

--------------------------------------------------------------------------------
FILE DESCRIPTION
--------------------------------------------------------------------------------

Below is an overview of the file and each function:

1) thermal_sizing_v1.m [MAIN SCRIPT]
   This is the project's entry point. It defines all input parameters (motor 
   specs, baseline geometry, CEA gas data, and structural limits). It 
   coordinates the workflow by calling geometry construction, TBC sizing, 
   and cooling channel optimization routines in sequence. 

2) function nozzle_geometry_section_array_plot
   Responsible for the spatial discretization (1D mesh) of the nozzle. Based 
   on the lengths and angles of the convergent and divergent sections, it 
   generates axial coordinates (x), local radius, area ratios, axial steps (dx), 
   and arc length elements (ds).

3) function cea_interpolations_plot
   A support function that takes thermodynamic data tables exported from 
   CEA (gamma, cp, Prandtl, viscosity, Temperature, Pressure, Mach) and 
   interpolates them onto the geometric node grid. It handles the 
   transition between equilibrium (convergent) and frozen (divergent) flow.

4) function mu_w_to_temperature
   Calculates the thermophysical properties of the cooling water (specifically 
   dynamic viscosity "mu" and Prandtl number "Pr") as a function of local 
   temperature. Used during iterative loops to update the water's Reynolds 
   and Nusselt numbers.

5) function tbc_thickness_calculation_plot
   Manages the thermal wall simulation. It first calculates the heat flux (q) 
   for a theoretical bare wall case to demonstrate the need for a TBC. It then 
   performs the sizing of the YSZ ceramic layer required to keep the copper 
   wall and bond coat temperatures within safe operational limits.

6) function water_channel_sizing_dp_min_plot
   Using the actual heat flux 
   (mitigated by the TBC), it calculates the optimal cooling channel geometry. 
   It explores a design space of possible configurations (height, width, number 
   of channels) to find the solution that minimizes total pressure drop (Delta P) 
   while satisfying geometric, velocity, and "No Boiling" constraints.

7) function water_channel_sizing_dp_min_SS
   A variant of the sizing function designed to export the entire solution 
   domain (Design Search Space). It is used to generate contour maps 
   (Delta P vs. Mass flow vs. Number of channels) to visually justify the 
   choice of the optimal configuration and highlight technical feasibility 
   boundaries.
