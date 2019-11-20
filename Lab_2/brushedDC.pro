// Brushed DC machine
// 4th year proyect BEAMS-ULB (2012-13): Finite Element and Experimental Analysis of DC machines
// Author: Olivier Candeur; Supervisors: J. Gyselinck, Y. Mollet
// GetDP-Gmsh model: R. Sabariego - November 2013, minor modifications October 2016

Include "brushedDC_data.geo";

DefineConstant[
  Flag_NL = { 1, Choices{0,1}, Name "Input/60Nonlinear BH-curve", ReadOnly 0, Visible 1}
];

Group {
  Stator_Core   = Region[{STATOR_CORE}]; // Ferromagnetic part of the stator
  Stator_Air    = Region[{STATOR_SLOT_OPENING}]; // Air in the stator
  Stator_Airgap = Region[{STATOR_AIRGAP}] // Airgap in the stator side
  Stator_Bnd_MB = Region[{STATOR_BND_MOVING_BAND}]; // Line limiting the airgap on the stator side, i.e. boundary of the moving band

  Rotor_Core   = Region[{ROTOR_CORE}]; // Ferromagnetic part of the rotor
  Rotor_Air    = Region[{****}]; // Air in the rotor
  Rotor_Airgap = Region[{ROTOR_AIRGAP}]; // Airgap in the rotor side
  Rotor_Bnd_MB = Region[{ROTOR_BND_MOVING_BAND}]; // Line limiting the airgap on the rotor side, i.e. boundary of the moving band

  MovingBand_PhysicalNb = #MOVING_BAND ; // Fictitious number for moving band, not in the geo file, this region is created and meshed by the FE solver

  Surf_bn0 = Region[{SURF_EXT}] ; // Lines where the boundary condition normal magnetic flux (bn) must be equal to zero

  Stator_CoilSides_P = Region[{STATOR_COILSIDES_P}] ; // Stator Coil regions with Positive current
  Stator_CoilSides_N = Region[{STATOR_COILSIDES_N}] ; // Stator coil regions with Negative current
  Stator_Winding = Region[ {Stator_CoilSides_P, Stator_CoilSides_N} ] ;

  // Rotor conductors
  Rotor_Slots_P = Region[{ROTOR_SLOTS_P}] ; // Rotor conductors with Positive current
  Rotor_Slots_N = Region[{ROTOR_SLOTS_N}] ; // Rotor conductors with Negative current
  Rotor_Winding = Region[ {Rotor_Slots_P, Rotor_Slots_N} ];

  Stator = Region[{ Stator_Core }] ;
  Rotor  = Region[{ Rotor_Core }] ;

  // Moving band:  with or without symmetry, the BND line of the rotor must be complete
  Stator_Bnd_MB = Region[{STATOR_BND_MOVING_BAND}]; // Inner boundary line of the air on the stator side
  Rotor_Bnd_MB  = Region[{ROTOR_BND_MOVING_BAND}];  // Outer boundary line of the airgap on the rotor side

  // Source regions
  Windings = Region[ {Stator_Winding, Rotor_Winding} ] ;

  // Movingband definition
  MB  = MovingBand2D[ MovingBand_PhysicalNb, Stator_Bnd_MB, Rotor_Bnd_MB, 1] ;

  // All air regions
  Air = Region[{ Rotor_Air, Rotor_Airgap, Stator_Air, Stator_Airgap, MB } ] ;

  // Complete domain
  Domain = Region[{ Air, Stator, Rotor, Stator_Winding, Rotor_Winding}];

  If(Flag_NL)
    DomainNL = Region[ {Stator_Core, Rotor_Core } ]; // Domain with nonlinear materials
    DomainL  = Region[ {Domain,-DomainNL} ]; // Domain with linear materials
  EndIf

}

// --------------------------------------------------------------------------
// --------------------------------------------------------------------------

Function {
  // Data for modeling a stranded inductor
  DefineConstant[
    Ie = { ***, Name "Input/51Ie stator field excitation current", Highlight "AliceBlue" },
    Ia = { ***, Name "Input/52Ia rotor armature current", Highlight "AliceBlue"}
  ] ;
  // Note that DefineConstant with Name definition is only used in the GUI, you could as well simply define
  // Ie = ***;
  // Ia = ***;

  Stator_CoilSide_Area[] = *** ; // Area of the cut of one of the stator windings
  Rotor_Slot_Area[] = *** ;      // Area of one of the rotor conductors

  nb_turns_stator_coils     = *** ; // number of turns per stator coil
  nb_conductors_rotor_slots = *** ; // number of conductors per rotor slot

  Ifac_stator_coils     = *** ; // fraction of Ie flowing in each turn of each field winding coil
  Ifac_rotor_conductors = *** ; // fraction of Ia flowing in each conductor of each rotor slot


  // Imposed current density in the windings
  // Define js[] for source regions in stator and rotor.
  // Use: js[#{***}] = Vector[***, ***,***];

}


// --------------------------------------------------------------------------
// --------------------------------------------------------------------------

Include "machine_magsta_a.pro" ;
