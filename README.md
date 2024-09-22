# Lumped-vortex panel method

## Introduction
This is the simplest model which can be found under the panel method category. Lumped (or _discrete_) vortex means that every panel has a discrete vortex located at the panel quarter chord point. This model is suited to analyze thin airfoils, since the mean curvature line is discretized in panels. The present code is designed in order to handle the multi-element problem.

This project is useful for students and curious of aerodynamics, since despite being very simple from the theory point of view, one can:
- Achieve strong insight on how to apply the fundamentals of aerodynamics.
- Understand how the panels and/or airfoils interact with each others.
- Take the first steps in programming (the code is written as much as possible as OOP).
- Investigate some numerical issues which arise when discretizing in a large number of panels.
- Have useful and estimated results for validation with other more sofisticated methods.

## Using the script
Once downloaded, you can try running `main.py` and the following will happen:
- Three windows showing:
  - Plot of two parabolic airfoils in tandem, called _Front_ and _Rear_ in the domain.
  - Plot of delta C<sub>p</sub> (difference between lower and upper C<sub>p</sub>) along relative chord (x/c) for _Front_ airfoil.
  - Plot of delta C<sub>p</sub> along relative chord for _Rear_ airfoil.
- Airfoils aerodynamic characteristics shown in terminal window:
  - Gamma: circulation associated to the specific airfoil.
  - L: lift (per unit span) generated by the specific airfoil.
  - C<sub>L</sub>: lift coefficient generated by the specific airfoil.
- A file called `results.dat` with the same results listed above.

### Setting an analysis
1. In line 7 of file `main.py` declare the incidence angle (in degrees) of the free stream velocity respect th _x_ axis, this angle is called _alpha_ as can be seen in the following image (from _Low-Speed Aerodynamics_ by Katz & Plotkin).
2. The freestream velocity magnitude is declared in line 8.
3. Airfoils are declared appending objects initialized with `Airfoil('name', 'coordinates_file.dat')` to `airfoils` array, as seen in lines 15 and 16 of file `main.py`.
4. Finally, run `main.py`.

![image](https://github.com/pzitelli84/discrete-vortex-panel-method/assets/8440605/74f7adc4-d11f-47ed-9bd6-cb4acae9e762)

### Defining airfoils
Each airfoil used in the problem must be defined separately in a file containing ordered nodes (rows), which columns represent the _x_ and _y_ coordinates. The files must be in the same directory as `main.py`. Take for instance `coord_parabolic_front_20.dat`, where you can see how coordinates are defined. The order of nodes must follow the orientarion shown in te above image, i.e.: each panel is defined with the frst node upstream and the second downstream.

When creating an object instance of `Airfoil` class, two parameters must be passed:
- A string which defines its name for identification purposes.
- A second string with the nodes coordinates file name.

Several coordinate files can be found inside `airfoils` directory, in their names the type of airfoil is described as well as the number of panels used. It is importan to remark that the nodes coordinates are measured in the global reference system, so when having a multi-element configuration, the location of every airfoil is taken from the nodes coordinates file.

### `main.py` overview
This is the skeleton of the script, where instructions to solve the problem are found. It uses methods defined in `utils.py`. The logics of the step-by-step can be followed straight forward and it is recommended not to modify the structure.

### `utils.py` overview
All the class definitions, functions and calculations are written in this file. The code is commented so the user can quickly understand what is being done by the script. This is fundamental when comparing the code to the math found in the theory.

## Theory background
As was mentioned in the introduction, this approach is the simplest one can find under the vast variety of panel methods flavours. Due to its simplicity, it is very easy to be coded. Since there is not any distribution of vorticity along the panel length but a discrete vortex attached to it, the formulation of the flow tangency boundary condition is found straight forward.

The induced velocities are found on every control point (located at the three-quarter panel length) pplying the Biot-Savart Law. Then, the dot product of the induced velocities times the unit normal vector of each panel is equated to zero. A linear algebraic equation system is obtained, where the vortex strengths of each panel are the unknowns.

Once the system is solved, the vortex strengths are used to calculate the C<sub>p</sub> difference along the airfoils, their individual lift using the Generalized Kutta-Joukowski Theorem as well as the individual C<sub>L</sub>.

### Numerical issues
An interesting issue with this type of methods is that, contrary to intuition, oscillations on the C<sub>p</sub> difference along the airfoils chord are observed when the airfoil is discretized in a larger number of panels. The reader is encouraged to investigate this phenomenon typically found in methods with discrete unknowns values instead of distributed ones. Katz & Plotkin give some insight on this.

### Recommended bibliography
For a broader and detailed explanation, the reader is recommended to study from:
- Low-Speed Aerodynamics, J. Katz & A. Plotkin.
- An Introduction to Theoretical and Computational Aerodynamics, J. Moran.
- Theoretical and Computational Aerodynamics, T. Sengupta.

## Future improvements
- Streamline plots.
