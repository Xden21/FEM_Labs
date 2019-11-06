/*
Tip: In Windows, you may use Notepad++ for editing this file, select C++ mode for seeing comments in color

Finite elements and electromagnetic fields (H9x34A), KU Leuven, R. V. Sabariego, Oct. 2017

2D FE modelling of an inductor with EI core, including
- linear and nonlinear (isotropic) material modelling
- static computation
- inductance computation
- force computation: Lorentz, finite-difference of co-energy, and virtual work (vw, virtual distortion of the FE mesh)

companion files:
- magsta_a.pro          : MVP finite element formulation
- BH.pro                : nonlinear material law

TO BE MODIFIED, this file:
- inductor2d.pro        : link between FE formulation and geometry: Groups, Functions

TO BE CREATED:
- inductor2d.geo        : generation of the geometry and the mesh
- inductor2d_data.geo   : parameters related to the geometry, mesh, material,...

==> this last file containing parameters is OPTIONAL, its creation is handy as it should include data common to the geo and the pro files.
If you prefer, you will have to repeat variables and values in different places, what is a source of errors

  For instance, let us define the numerical variable (note that it is not a string):
   CORE=1000;
   will appear in the geo-file as
   Physical Surface(CORE) ={***};
   and then in the pro-file as
   Core = Region[{CORE}]; // For the language interpreter, "Core" is a string, not CORE.

   Note that Gmsh allows you to define the string, and that helps checking the physical regions:
   Physical Surface("Core", CORE) ={***};

*/

// Uncomment the following line if you have created the data OPTIONAL file
// Include "inductor2d_data.geo";

psim='Simulation param./'; // For creating a tree in the GUI (not relevant for the computations)

DefineConstant[
  Flag_BC_Type = {1, Choices{0='Neumann',1='Dirichlet'},
    Name StrCat[psim,'20Boundary condition at infinity']}
  Flag_NL = { 0, Choices{0,1},Name StrCat[psim,'40Nonlinear BH-curve']}
  Flag_Infinity = {0, Choices{0,1},
    Name StrCat[psim,'21Use shell transformation to infinity']}
];


Group {
  // Link with the Physical Regions in geo-file
  // Note that the following regions are related to the physical properties
  // and could comprise more than one geometrical entity
  Core = Region[ {***} ];

  Inds   = Region[ {***} ];

  Air    = Region[ {***} ];
  AirInf = Region[ {***} ];


 // Surfaces for imposing boundary conditions
 If(Flag_BC_Type==1)
   Surf_Inf = Region[ {***} ]; // for imposing BC az=0 or not (Dirichlet vs Neumann)
 EndIf

 Surf_bn0 = Region[ {***} ]; // when modelling only half of the geometry (with az=0 on the symmetry axis)

 SkinDomain_Moving = Region[{***}] ; // for force computation on Icore using the VW method
 NumCores = 1; // Only ICore consider for force computation
 SkinCore_1 = Region[{SkinDomain_Moving}];

 }

Function {
  Freq = 50. ;

  SymmetryFactor = ***;
  LengthAlongZ = ***;

  // Consider a stranded inductor
  Irms =  *** ; // Current (rms) [A]
  II = *** ; // Current magnitude

  NbWires[]  = ***; // Number of turns
  SurfCoil[] = ***; // Surface area of the coil
  Idir[]     = ***; // direction of imposed current: +1 coming out of plane, -1 going in the plane

  vDir[]   = Vector[ 0, 0, Idir[]] ; // vector indicating the direction

  js0[] = *** ; //imposed current density

  // Material properties
  mu0 = 4.e-7 * Pi ;

  // Use if linear problemin
  sigma_coil = *** ;// Conductivity of copper [S/m]
  mur_fe =  *** ; //Core relative permeability, only use in case of linear problem

  // Cylindrical shell radii
  Val_Rint = ***; // inner radius [m]
  Val_Rext = ***; // outer radius [m]

 }

 Include 'magsta_a.pro' ;
