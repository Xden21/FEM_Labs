// the following general group names are defined,
//      and may then be formally used in the rest of the pro-file without prior assignment

Group {
   DefineGroup[
     DomainB,   // domain with stranded coils (allowing for circuit coupling), 'B' of 'bobiné'
     DomainS,   // domain with conductors in which current density is imposed directly (not used here), 'S' of 'source'
     DomainInf, // domain where transformation to infinity takes places
     DomainL,   // domain with linear material (field-independent permeability), isotropic only here
     DomainNL,  // domain with nonlinear materials (field-dependent permeability), isotropic only here
     Domain,    // complete domain

     Surf_bn0, // imposing normal component of induction, equivalent to impose az value, e.g. accounting for symmetry
     Surf_Inf  //for imposing BC az = 0 or not (Dirichlet vs. Neumann)
  ];
}

Function{
  DefineConstant[ //Predefined values, that will be changed if defined somewhere else in the pro-file
    Nb_max_iter        = 30,
    relaxation_factor  = 1,
    stop_criterion     = 1e-5,
    // FillFactor_Winding = {1, Name '{Simulation param./4Coil/3Fill factor', Highlight 'AliceBlue', Visible 1},
    // Factor_R_3DEffects = {1, Name '{Simulation param./4Coil/4 3D factor', Highlight 'AliceBlue', Visible 1},
    // Increasing the resistance by 50% == 1.5
    II, VV,
    Flag_NL_law_Type = 0,
    po    = '{Output/0'
    poL   = '{Output/11Inductance/1'
    poLA  = '{Output/10Inductance Analytical/1'
    poF   = '{Output/21Force/2'
    poFA  = '{Output/20Force Analytical/2'
  ];

}

// Typical file extensions
ExtGmsh     = '.pos'; // Default extension for Gmsh files
ExtGnuplot  = '.dat'; // Just a text file, other people prefer .txt :-)

Include 'BH.pro'; // nonlinear BH characteristic of magnetic material

Group {
  DomainS   = Region[ {Inds} ];

  If(Flag_Infinity)
    DomainInf = Region[ {AirInf} ];
  EndIf

  Domain = Region[ {Air, AirInf, Inds, Core} ];

  If(!Flag_NL)
    DomainNL = Region[ {} ];
    DomainL  = Region[ {Domain} ];
  EndIf
  If(Flag_NL)
    DomainNL = Region[ {Core} ];
    DomainL  = Region[ {Domain, -DomainNL} ];
  EndIf

  DomainDummy = #12345 ; // Dummy region number for postpro with functions
}

Function {
  nu [#{Air, AirInf, Inds}]  = 1./mu0 ;

  If(!Flag_NL)
    nu [#{Core}]  = 1/(mur_fe*mu0) ;
    DefineFunction[ dhdb_NL ];
  EndIf
  If(Flag_NL)
    nu [ DomainNL ] = nu_EIcore[$1] ;
    dhdb_NL [ DomainNL ] = dhdb_EIcore_NL[$1];
    // nu [ DomainNL ] = nu_brauer[$1] ; // Analytical law approximated by Matlab
    // dhdb_NL [ DomainNL ] = dhdb_brauer_NL[$1];
  EndIf

  sigma[#{Inds}] = sigma_coil ;
  rho[] = 1/sigma[] ;

  // energy and co-energy density (J/m^3), functions of the real b-vector
  w_mag[DomainL]    = 0.5*nu[]*SquNorm[$1] ;  // linear materials
  w_mag_co[DomainL] = w_mag[$1] ;
  w_mag[DomainNL]   = 0.5*k1*SquNorm[$1]+0.5*k2/k3*(Exp[k3*SquNorm[$1]]-1) ; // nonlinear materials
  w_mag_co[DomainNL]= nu[$1]*SquNorm[$1] - w_mag[$1] ;

  // Maxwell tensor
  Tmax[] = (SquDyadicProduct[$1]-SquNorm[$1]*TensorDiag[0.5,0.5,0.5])/mu0 ;
}

//-------------------------------------------------------------------------------------

Jacobian {
  { Name Vol;
    Case {
      { Region DomainInf ; Jacobian VolSphShell{Val_Rint, Val_Rext} ; }
      { Region All ; Jacobian Vol; }
    }
  }
  { Name Sur;
    Case {
      { Region All ; Jacobian Sur; }
    }
  }
}

Integration {
  { Name I1 ; Case {
      { Type Gauss ;
        Case {
          { GeoElement Triangle   ; NumberOfPoints  6 ; }
	  { GeoElement Quadrangle ; NumberOfPoints  4 ; }
	  { GeoElement Line       ; NumberOfPoints  13 ; }
        }
      }
    }
  }

  { Name I_1p ; Case { // Use in post-processing for integrating quantities that are constant per element
      { Type Gauss ;
        Case {
          { GeoElement Line       ; NumberOfPoints  1 ; }
          { GeoElement Triangle   ; NumberOfPoints  1 ; }

        }
      }
    }
  }
}

//-------------------------------------------------------------------------------------

Constraint {

  { Name MVP_2D ;
    Case {
      { Region Surf_Inf ; Type Assign ; Value 0. ; }
      { Region Surf_bn0 ; Type Assign ; Value 0. ; }
    }
  }

  For i In {1:NumCores}
    { Name auxDir~{i} ;
      Case {
        { Region SkinCore~{i} ; Value -1. ; }
      }
    }
  EndFor

}

//-----------------------------------------------------------------------------------------------

FunctionSpace {
  // Magnetic Vector Potential
  { Name Hcurl_a_2D ; Type Form1P ;
    BasisFunction {
      { Name se1 ; NameOfCoef ae1 ; Function BF_PerpendicularEdge ;
        Support Region[{Domain, SkinDomain_Moving}] ; Entity NodesOf [ All ] ; }
   }
    Constraint {
      { NameOfCoef ae1 ; EntityType NodesOf ; NameOfConstraint MVP_2D ; }
    }
  }

  // auxiliary field on layer of elements touching region of interest, for the
  // accurate integration of the Maxwell stress tensor (using the gradient of
  // this field)
  For i In {1:NumCores}
    { Name auxDir~{i} ; Type Form0 ;
      BasisFunction {
        { Name sn ; NameOfCoef un ; Function BF_GroupOfNodes ;
          Support Air ; Entity GroupsOfNodesOf[ SkinCore~{i} ] ; }
      }
      Constraint {
        { NameOfCoef un ; EntityType GroupsOfNodesOf ; NameOfConstraint auxDir~{i} ; }
      }
    }
  EndFor

}

//-----------------------------------------------------------------------------------------------

Formulation {

  { Name MagSta_a_2D ; Type FemEquation ;
    Quantity {
      { Name a  ; Type Local  ; NameOfSpace Hcurl_a_2D ; }

      For i In {1:NumCores}
        { Name un~{i} ; Type Local ; NameOfSpace auxDir~{i} ; }
      EndFor
    }

    Equation {
      Galerkin { [ nu[{d a}] * Dof{d a}  , {d a} ] ;
        In Domain ; Jacobian Vol ; Integration I1 ; }
      Galerkin { JacNL [ dhdb_NL[{d a}] * Dof{d a} , {d a} ] ;
        In DomainNL ; Jacobian Vol ; Integration I1 ; }

      Galerkin { [ -js0[] , {a} ] ;
        In DomainS ; Jacobian Vol ; Integration I1 ; }

      // dummy term to define dofs for fully fixed space
      For i In {1:NumCores} // dummy term to define dofs for fully fixed space
        Galerkin { [ 0 * Dof{un~{i}} , {un~{i}} ] ;
          In Domain ; Jacobian Vol ; Integration I1 ; }
      EndFor
    }
  }

}

//-----------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------

Resolution {

  { Name Analysis ;
    System {
      { Name A ; NameOfFormulation MagSta_a_2D ; }
    }
    Operation {

      If(!Flag_NL)
        Generate[A] ; Solve[A] ;
      EndIf
      If(Flag_NL)
        IterativeLoop[Nb_max_iter, stop_criterion, relaxation_factor]{
          GenerateJac[A] ; SolveJac[A] ; }
      EndIf
      SaveSolution[A] ;

      PostOperation[Get_LocalFields] ;
      PostOperation[Get_GlobalQuantities] ;
    }
  }
}

//-----------------------------------------------------------------------------------------------

PostProcessing {

  { Name MagSta_a_2D ; NameOfFormulation MagSta_a_2D ;
    PostQuantity {
      { Name a  ; Value { Term { [ {a} ] ; In Domain ; Jacobian Vol ; } } }
      { Name az ; Value { Term { [ CompZ[{a}] ] ; In Domain ; Jacobian Vol ; } } }

      { Name b  ; Value { Term { [ {d a} ] ; In Domain ; Jacobian Vol ; } } }
      { Name nb  ; Value { Term { [ Norm[{d a}] ] ; In Domain ; Jacobian Vol ; } } }

      { Name h ; Value { Term { [ nu[{d a}] * {d a} ] ; In Domain ; Jacobian Vol ; } } }

      { Name js0 ; Value { Term { [ js0[] ] ; In DomainS ; Jacobian Vol ; } } }

      // magnetic energy and co-energy density, J/m^3
      { Name w_mag ;    Value { Term { [ w_mag[{d a}]    ] ; In Domain ; Jacobian Vol ; } } }
      { Name w_mag_co ; Value { Term { [ w_mag_co[{d a}] ] ; In Domain ; Jacobian Vol ; } } }
      // Magnetic energy and co-energy
      { Name W_mag ; Value {
          Integral { [ SymmetryFactor*LengthAlongZ * w_mag[{d a}] ] ;
            In Domain ; Jacobian Vol ; Integration I1 ; } } }
      { Name W_mag_co ; Value {
          Integral { [ SymmetryFactor*LengthAlongZ * w_mag_co[{d a}] ] ;
            In Domain ; Jacobian Vol ; Integration I1 ; } } }

      { Name Flux ; Value { // Flux linkage computation
          Integral { [ SymmetryFactor*LengthAlongZ*Idir[]*NbWires[]/SurfCoil[]* CompZ[{a}] ] ;
            In Inds  ; Jacobian Vol ; Integration I1 ; }
        }
      }

      // Inductance computation from previously computed global quantities
      // Flux stored in register #11 => variable $flux
      // Magnetic energy, W_mag, stored in register #22 => variable $wmag
      // Magnetic co-energy, W_mag_co, stored in register #23 => variable $wmag_co

      { Name L_from_Flux ; Value { Term { Type Global; [ $flux*1e3/II ] ; In DomainDummy ; } } }
      { Name L_from_W_mag ; Value { Term { Type Global; [ 2*$wmag*1e3/(II*II) ] ; In DomainDummy ; } } }
      { Name L_from_W_mag_co ; Value { Term { Type Global; [ 2*$wmag_co*1e3/(II*II) ] ; In DomainDummy ; } } }

      // Force (N) computation by Maxwell stress tensor (exterior normal automatic from geo, surface orientation matters)
      { Name Fskin  ; Value { Integral { [ Lz * SymmetryFactor * Tmax[{d a}] * Normal[] ] ;
            In SkinDomain_Moving ; Jacobian Sur ; Integration I_1p ; } } } // Force_skin_mst stored in register #54
      { Name Fskinx ; Value { Term { Type Global; [ CompX[$fskin] ] ; In DomainDummy ; } } }
      { Name Fskiny ; Value { Term { Type Global; [ CompY[$fskin] ] ; In DomainDummy ; } } }

      // Force computation using MST
      For i In {1:NumCores}
        // Aux quantity for computing correctly the normal
        { Name un~{i} ; Value { Local { [ {un~{i}} ] ; In Domain ; Jacobian Vol ; } } }
        { Name dun~{i} ; Value { Local { [ {d un~{i}} ] ; In Domain ; Jacobian Vol ; } } }

        { Name f~{i}  ; Value { Integral { [ LengthAlongZ *SymmetryFactor * Tmax[{d a}] * {d un~{i}} ] ;
              In Air ; Jacobian Vol ; Integration I1 ; } } }
        { Name fx~{i}  ; Value { Integral { [ CompX[LengthAlongZ *SymmetryFactor * Tmax[{d a}] * {d un~{i}}] ] ;
              In Air ; Jacobian Vol ; Integration I1 ; } } }
        { Name fy~{i}  ; Value { Integral { [ CompY[LengthAlongZ *SymmetryFactor * Tmax[{d a}] * {d un~{i}}] ] ;
              In Air ; Jacobian Vol ; Integration I1 ; } } }
      EndFor

      // checking Normal definion
      { Name N_line  ; Value { Term { [ Normal[] ] ; In SkinDomain_Moving ; Jacobian Vol ; } } }


    }//PostQuantity
  }// MagStaDyn_a_2D
}// PostProcessing


//-----------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------

PostOperation Get_LocalFields UsingPost MagSta_a_2D {
  Print[ un_1, OnElementsOf Domain, File StrCat['un',ExtGmsh], LastTimeStepOnly ];
  Print[ dun_1, OnElementsOf Domain, File StrCat['dun',ExtGmsh], LastTimeStepOnly ];

  Print[ b,  OnElementsOf Domain, File StrCat['b',ExtGmsh], LastTimeStepOnly ] ;
  Print[ nb,  OnElementsOf Domain, File StrCat['nb',ExtGmsh], LastTimeStepOnly ] ;
  Print[ az, OnElementsOf Domain, File StrCat['a',ExtGmsh], LastTimeStepOnly ];

  Echo[Str['View[PostProcessing.NbViews-1].RangeType = 3;' ,// per timestep
      'View[PostProcessing.NbViews-1].NbIso = 25;',
      'View[PostProcessing.NbViews-1].IntervalsType = 1;' // iso values
    ], File 'az.opt'];
}

PostOperation Get_GlobalQuantities UsingPost MagSta_a_2D {
  Print[ Flux[Inds], OnGlobal, Format TimeTable,
    File > StrCat['Flux',ExtGnuplot], LastTimeStepOnly, StoreInVariable $flux,
    SendToServer StrCat[po,'40Flux [Wb]'],  Color 'LightYellow' ]; //Wheat

  Print[ W_mag[Domain], OnGlobal, Format TimeTable,
    File > StrCat['W_mag',ExtGnuplot], LastTimeStepOnly, StoreInVariable $wmag,
    SendToServer StrCat[po,'41Mag. energy [J]'],  Color 'LightYellow' ];
  Print[ W_mag_co[Domain], OnGlobal, Format TimeTable,
    File > StrCat['W_mag_co',ExtGnuplot], LastTimeStepOnly, StoreInVariable $wmag_co,
    SendToServer StrCat[po,'42Mag. co-energy [J]'],  Color 'LightYellow' ];

  Print[ L_from_Flux, OnRegion DomainDummy, Format Table, LastTimeStepOnly, File StrCat['L_flux',ExtGnuplot],
    SendToServer StrCat[poL,'50L from flux [mH]'], Color 'LightYellow' ];
  Print[ L_from_W_mag, OnRegion DomainDummy, Format Table, LastTimeStepOnly, File StrCat['L_energy',ExtGnuplot],
    SendToServer StrCat[poL,'51L from mag. energy [mH]'], Color 'LightYellow' ];
  Print[ L_from_W_mag_co, OnRegion DomainDummy, Format Table, LastTimeStepOnly, File StrCat['L_coenergy',ExtGnuplot],
    SendToServer StrCat[poL,'52L from mag. co-energy [mH]'], Color 'LightYellow' ];

  // Maxwell stress tensor on a line
  Print[ Fskin[SkinDomain_Moving], OnGlobal, Format TimeTable,
    File > StrCat['Fmst_skin',ExtGnuplot], StoreInVariable $fskin, LastTimeStepOnly ];
  Print[ Fskinx, OnRegion DomainDummy, Format Table, LastTimeStepOnly,
    File 'tmp.dat', SendToServer StrCat[poF,'61Fx (line) [N]'], Color 'LightYellow'];
  Print[ Fskiny, OnRegion DomainDummy, Format Table, LastTimeStepOnly,
    File 'tmp.dat', SendToServer StrCat[poF,'62Fy (line) [N]'], Color 'NavajoWhite'];

  // Maxwell stress tensor on a surface, exact computation of normal
  For i In {1:NumCores}
    //Print[ un~{i}, OnElementsOf Domain, File Sprintf('un_%g.pos',i) ];
    Print[ f~{i}[Air],  OnGlobal, Format TimeTable,
           File > StrCat[Sprintf('Fmst%g',i),ExtGnuplot], StoreInVariable $fmax, LastTimeStepOnly ];
    Print[ fx~{i}[Air], OnGlobal, Format Table, LastTimeStepOnly,
      File 'tmp.dat', SendToServer StrCat[poF,'65Fx (surf) [N]'], Color 'Ivory'  ];
    Print[ fy~{i}[Air], OnGlobal, Format Table,
      File 'tmp.dat', SendToServer StrCat[poF,'66Fy (surf) [N]'], Color 'NavajoWhite'  ];
  EndFor
 }


DefineConstant[
   // always solve this resolution (it's the only one provided)
  R_ = {'Analysis', Name 'GetDP/1ResolutionChoices', Visible 0}, // Nothing to worry about
   // set some command-line options for getdp
  C_ = {'-solve -v 3 -v2', Name 'GetDP/9ComputeCommand', Visible 1},
  // don't do the post-processing pass (some already included in Resolution)
  P_ = {'', Name 'GetDP/2PostOperationChoices', Visible 1} // nothing to worry about ...
];
