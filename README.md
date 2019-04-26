# rotor-split
Computes frequencies (rad/s or Hz), modeshapes and NAVMI factor of free-clamped rings &amp; clamped circular plates rotating
in a fluid for a given number of nodal diameters and circles; the geometry is customizable. Also allows for simple Dirac forced
response analysis.

%% Parameters\\
Input geometry and material properties as well as considered mode

%% k calculation
Computes the natural frequency of the structure in vacuum, according to the model from Leissa (1969)
The modeshape is assumed to be a sum of Bessel functions and applying the boundary conditions gives a system of equation
Solving the det=0 gives k

%% Omega calculation
Natural frequency in vacuum is obtained from k

%% Modeshape
Solving the boundary condition equations (with an arbitrary norm) gives the coefficient for the modeshape

%% Displacement W
Plots the modeshape

%% Added mass M
Computes H, which is derived from the model of Amabili et al. (1996) (using the Hankel transform method). It is a function
of the modeshape and consists in non trivial integrals that need to be solved numerically
The added mass can then be determined from this expression

%% NAVMI factor
Computes the reference kinetic energy of the disk and of the fluid (both above and below the disk)
Computes the AVMI factor as the ratio of these kinetic energies
Computes the natural frequencies of the rotating disk in fluid as a function of the rotation speed, the AVMI factor and the
natural frequency in vacuum (cf. analytical model from Max Louyot)

%% OmegaD range plot
Plots natural frequencies for a range of rotation speeds

%% Forced response -- Dirac excitation
Plots the forced response amplitude of a point on the disk when excited with a Dirac
