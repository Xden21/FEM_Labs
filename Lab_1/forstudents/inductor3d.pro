Include 'inductor3d_data.geo';

Dir='res/';
ExtGmsh     = '.pos';
ExtGnuplot  = '.dat';

TREE_COTREE_GAUGE=0;
COULOMB_GAUGE=1;

DIVJ0_NONE = 0;
DIVJ0_WEAK = 1;

DefineConstant[
  Flag_BC_Type = {1, Choices{0='Neumann',1='Dirichlet'}, ReadOnly (Flag_Infinity==1),
    Name StrCat[psim,'20Boundary condition at infinity'], Highlight 'Blue', Visible 0},

  Flag_NL = { 0, Choices{0,1},
    Name StrCat[psim,'40Nonlinear BH-curve'], Visible 1}
  Flag_NL_brauer = { 0, Choices{0,1},
    Name StrCat[psim,'41Analytical Brauer-type law?'], Visible Flag_NL, Highlight 'Blue'}

  Flag_GaugeType = { 0, Choices{0='Tree-cotree gauge', 1='Coulomb gauge'},
    Name StrCat[psim,'42Type of gauge'], Highlight 'Blue', Visible 1 }

  Flag_DivJ_Zero = { DIVJ0_NONE, Choices{ DIVJ0_NONE = 'None', DIVJ0_WEAK = 'Weak'},
    Name StrCat[psim,'30Constraint div j = 0'],
    Help Str['None: direct interpolation of js0[]',
      'Weak: Use scalar potential xis for weakly ensuring div j = 0.',
      'Strong: Use Hcurl source field hs with curl hs = j, for div j = 0;'],
    Highlight 'Blue'}
];

Group {
  Core = Region[ {ECORE, ICORE} ];
  Bnd_Icore = Region[{SKINICORE}];
  Bnd_Ecore = Region[{SKINECORE}];
  SkinDomain_Moving = Region[ {SKINICORE} ];
  Inds  = Region[ {COIL} ];
  SkinInds = Region[ {SKINCOIL} ];

  Air  = Region[ {AIR, AIRGAP} ];
  AirInf = Region[ AIRINF ];

 // Surfaces for imposing boundary conditions
 If(Flag_BC_Type==1)
   Surf_Inf = Region[ SURF_AIROUT ];
 EndIf

 Surf_bn0 = Region[ {CUT_YZ, CUT_XY} ];

 NumRegionsForceComputation = 1;
 SkinRegionForceComputation_1 = Region[ {SKINICORE} ];
}

Function {
  Freq = 50. ;

  DefineConstant[
    Irms = { IA, Min 1, Max 2*IA, Step 1,
      Name StrCat[psim,'4Coil/0Current (rms) [A]'], Highlight 'AliceBlue'}
      NbWires = { Nw, Name StrCat[psim,'4Coil/1Number of turns'],
      Highlight 'AliceBlue'}
  ];

  II = Irms *Sqrt[2] ;
  NbWires[]  = NbWires;
  SurfCoil[] = wcoil*hcoil;
  vDir[] = -( // change of sign for coherence with 2D model
    (X[] >  wcoreE && Z[] <  -Lz/2) ? Vector[ Sin[ Atan2[Z[]+Lz/2, X[]-wcoreE]#1], 0, -Cos[#1]]:
    (X[] >= wcoreE && Z[] >= -Lz/2) ? Vector[ 0, 0, -1] : Vector[ -1, 0, 0] );

  js0[] = II*NbWires[]/SurfCoil[] * vDir[] ; //DomainS


  // Material properties
  mu0 = 4.e-7 * Pi ;

  DefineConstant[
    sigma_coil = { sigma_cu, Name StrCat[psim,'4Coil/5Conductivity'], Units 'S/m', Highlight 'AliceBlue'},
    mur_fe = { 1000., Min 100, Max 1000, Step 100,
      Name StrCat[psim,'42Core relative permeability'], Highlight 'AliceBlue',
      Visible (!Flag_NL)} // used if linear material characteristic
  ];

}

Include 'magsta_a_js0_3d.pro' ;
