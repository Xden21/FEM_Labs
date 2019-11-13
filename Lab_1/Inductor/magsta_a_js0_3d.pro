// --------------------------------------------------------------------------
// --------------------------------------------------------------------------

Group {
  DefineGroup[
    Domain, DomainL, DomainNL, DomainS, DomainInf,
    SkinDomain_Moving, Surf_bn0, Surf_Inf, Surf_FixedMVP
  ] ;
}


Function {
  DefineFunction[
    mu, nu, sigma, rho, js, dhdb_NL
  ] ;

  DefineConstant[
    Nb_max_iter        = 30,
    relaxation_factor  = 1,
    stop_criterion     = 1e-5,
    po       = "{Output/0"
    poL      = "{Output/1Inductance/"
    poF      = "{Output/2Force/"
  ];

}
Include "BH.pro"; // nonlinear BH caracteristic of magnetic material

Group {

  DomainS   = Region[ {Inds} ];
  SkinDomainS = Region[ {SkinInds} ];

  Domain = Region[ {Air, AirInf, Core, Inds} ];

  If(Flag_Infinity)
    DomainInf = Region[ {AirInf} ];
  EndIf

  If(!Flag_NL)
    DomainNL = Region[ {} ];
    DomainL  = Region[ {Domain} ];
  EndIf
  If(Flag_NL)
    DomainNL = Region[ {Core} ];
    DomainL  = Region[ {Domain,-DomainNL} ];
  EndIf

  DomainDummy = #12345 ; // Dummy region number for postpro with functions

  Surf_FixedMVP = Region[{ Surf_bn0, Surf_Inf }];
}

Function {
  nu [#{Air, AirInf, Inds}]  = 1./mu0 ;

  If(!Flag_NL)
    nu [#{Core}]  = 1/(mur_fe*mu0) ;
  EndIf
  If(Flag_NL)
    If(!Flag_NL_brauer)
      nu [ DomainNL ] = nu_EIcore[$1] ;
      dhdb_NL [ DomainNL ] = dhdb_EIcore_NL[$1];
    EndIf
    If(Flag_NL_brauer)
      nu [ DomainNL ] = nu_EIcore_a[$1] ;
      dhdb_NL [ DomainNL ] = dhdb_EIcore_a_NL[$1];
    EndIf
  EndIf

  sigma[#{Inds}] = sigma_coil ;
  rho[] = 1/sigma[] ;

   // energy and co-energy density (J/m^3), functions of the real b-vector
  w_mag[DomainL]    = 0.5*nu[]*SquNorm[$1] ;  // linear materials
  w_mag_co[DomainL] = w_mag[$1] ;
  w_mag[DomainNL]   = 0.5*k1*SquNorm[$1]+0.5*k2/k3*(Exp[k3*SquNorm[$1]]-1) ; // nonlinear materials
  w_mag_co[DomainNL]= nu[$1]*SquNorm[$1] - w_mag[$1] ;

  // Maxwell stress tensor => used in force computation
  Tmax[] = (SquDyadicProduct[$1]-SquNorm[$1]*TensorDiag[0.5,0.5,0.5])/mu0 ;

}


// --------------------------------------------------------------------------

Jacobian {
  { Name Vol ;
    Case { { Region DomainInf ; Jacobian VolSphShell {Val_Rint, Val_Rext} ; }
           { Region All ;       Jacobian Vol ; }
    }
  }
  { Name Sur ;
    Case { { Region All ; Jacobian Sur ; }
    }
  }
}

Integration {
  { Name II ;
    Case {
      {
	Type Gauss ;
	Case {
	  { GeoElement Triangle    ; NumberOfPoints  4 ; }
	  { GeoElement Quadrangle  ; NumberOfPoints  4 ; }
	  { GeoElement Tetrahedron ; NumberOfPoints  4 ; }
	  { GeoElement Hexahedron  ; NumberOfPoints  6 ; }
	  { GeoElement Prism       ; NumberOfPoints  21 ; }
	  { GeoElement Line        ; NumberOfPoints  4 ; }
	}
      }
    }
  }

   { Name I_1p ; Case {
      { Type Gauss ;
	Case {
	  { GeoElement Triangle    ; NumberOfPoints  1 ; }
	  { GeoElement Tetrahedron ; NumberOfPoints  1 ; }
	}
      }
    }
  }

}

// --------------------------------------------------------------------------

Constraint {
  // av - formulation
  { Name MVP_3D ;
    Case {
      { Region Surf_bn0 ; Type Assign ; Value 0. ; }
      { Region Surf_Inf ; Type Assign ; Value 0. ; }
    }
  }

  // Constraint on the auxiliary basis function for the correct computation of the force using the MST
  // Negative sign refers to the exterior normal at the skin of the considered region
  For i In {1:NumRegionsForceComputation}
    { Name auxDir~{i} ;
      Case {
        { Region SkinRegionForceComputation~{i} ; Value -1. ; }
      }
    }
  EndFor
}

Group {
  Surf_a_NoGauge = Region [ {Surf_FixedMVP} ] ;
}

Constraint {

  { Name GaugeCondition_a ; Type Assign ;
    Case {
      If(Flag_GaugeType==TREE_COTREE_GAUGE)
        { Region Domain ; SubRegion Surf_a_NoGauge ; Value 0. ; }
      EndIf
      }
  }

  { Name xi_fixed ; Type Assign ;
    Case {
      { Region Surf_FixedMVP ; Value 0. ; }
    }
  }


}

FunctionSpace {

  // Magnetic vector potential a (b = curl a)
  { Name Hcurl_a_3D ; Type Form1 ;
    BasisFunction {// a = a_e * s_e
      { Name se ; NameOfCoef ae ; Function BF_Edge ;
        Support Region[{Domain, SkinDomain_Moving}] ; Entity EdgesOf[ All ] ; }
    }
    Constraint {
      { NameOfCoef ae  ; EntityType EdgesOf ; NameOfConstraint MVP_3D ; }
      { NameOfCoef ae  ; EntityType EdgesOfTreeIn ; EntitySubType StartingOn ;
        NameOfConstraint GaugeCondition_a ; }
    }
  }

  // auxiliary field on layer of elements touching region of interest, for the
  // accurate integration of the Maxwell stress tensor (using the gradient of
  // this field)
  For i In {1:NumRegionsForceComputation}
    { Name auxDir~{i} ; Type Form0 ;
      BasisFunction {
        { Name sn ; NameOfCoef un ; Function BF_GroupOfNodes ;
          Support Air ; Entity GroupsOfNodesOf[ SkinRegionForceComputation~{i} ] ; }
      }
      Constraint {
        { NameOfCoef un ; EntityType GroupsOfNodesOf ; NameOfConstraint auxDir~{i} ; }
      }
    }
  EndFor


  // scalar potential for Coulomb gauge: a orthogonal to grad(Xi)
  { Name H_xi ; Type Form0 ;
    BasisFunction {
      { Name sn ; NameOfCoef an ; Function BF_Node ;
        Support Region[{Domain}] ; Entity NodesOf[ All ] ; }
    }
    Constraint {
      { NameOfCoef an ; EntityType NodesOf ; NameOfConstraint xi_fixed ; }
    }
  }

  // correcting source interpolation js0[] so that (weakly) div j = 0
  { Name H_xi_divj0 ; Type Form0 ;
    BasisFunction {
      { Name sn ; NameOfCoef an ; Function BF_Node ;
        Support Region[{DomainS, SkinDomainS}] ; Entity NodesOf[ All ] ; }
    }
    Constraint {
      { NameOfCoef an ; EntityType NodesOf ; NameOfConstraint xi_fixed ; }
    }
  }


}

//---------------------------------------------------------------------------------------------

Formulation {

  { Name DivJ0 ; Type FemEquation ;
    Quantity {
      { Name xi; Type Local ; NameOfSpace H_xi_divj0 ; }
    }
    Equation {
      Galerkin { [ js0[] , {d xi} ] ;
        In Domain ; Jacobian Vol ; Integration II ; }
      Galerkin { [ -Dof{d xi} , {d xi} ] ;
        In Domain ; Jacobian Vol ; Integration II ; }
    }
  }

  { Name MagSta_a_js0_3D ; Type FemEquation ;
    Quantity {
      { Name a  ; Type Local ; NameOfSpace Hcurl_a_3D ; }
      { Name xi ; Type Local ; NameOfSpace H_xi ; } // Coulomb gauge
      { Name xis ; Type Local ; NameOfSpace H_xi_divj0 ; } // div j=0

      For i In {1:NumRegionsForceComputation}
        { Name un~{i} ; Type Local ; NameOfSpace auxDir~{i} ; }
      EndFor
    }

    Equation {
      Galerkin { [ nu[{d a}] * Dof{d a} , {d a} ] ;
        In Domain ; Jacobian Vol ; Integration II ; }
      Galerkin { JacNL [ dhdb_NL[{d a}] * Dof{d a} , {d a} ] ;
        In DomainNL ; Jacobian Vol ; Integration II ; }

      Galerkin { [ -js0[], {a} ] ;
        In  DomainS ; Jacobian Vol ; Integration II ; }


      If(Flag_DivJ_Zero == DIVJ0_WEAK)
        Galerkin { [ {d xis}, {a} ] ;
          In Domain ; Jacobian Vol ; Integration II ; }
      EndIf

      If(Flag_GaugeType==COULOMB_GAUGE)
        Galerkin { [ Dof{a}, {d xi} ] ;
          In  Domain ; Jacobian Vol ; Integration II ; }
        Galerkin { [ Dof{d xi}, {a} ] ;
          In  Domain ; Jacobian Vol ; Integration II ; }
      EndIf

      // dummy term to define dofs for fully fixed space
      For i In {1:NumRegionsForceComputation}
        Galerkin { [ 0 * Dof{un~{i}} , {un~{i}} ] ;
          In Domain ; Jacobian Vol ; Integration II ; }
      EndFor
    }
  }

}

Resolution {

  { Name Analysis ;
    System {
      { Name Sys ; NameOfFormulation MagSta_a_js0_3D ; }
      If(Flag_DivJ_Zero == DIVJ0_WEAK)
        { Name Sys_DivJ0 ; NameOfFormulation DivJ0 ; }
      EndIf
    }
    Operation {
      CreateDir["res/"] ;

      If(Flag_DivJ_Zero == DIVJ0_WEAK)
        Generate[Sys_DivJ0]; Solve[Sys_DivJ0]; SaveSolution[Sys_DivJ0];
      EndIf

      If(!Flag_NL)
        Generate[Sys]; Solve[Sys];
      EndIf
      If(Flag_NL)
        IterativeLoop[Nb_max_iter, stop_criterion, relaxation_factor]{
          GenerateJac[Sys]; SolveJac[Sys]; }
      EndIf

      SaveSolution[Sys];

      PostOperation[Get_LocalFields];
      PostOperation[Get_GlobalQuantities];
    }
  }

}

//-----------------------------------------------------------------------------------------------

PostProcessing {

  { Name MagSta_a_3D ; NameOfFormulation MagSta_a_js0_3D ;
    PostQuantity {
      { Name a ; Value { Term { [ {a} ]          ; In Domain ; Jacobian Vol ; } } }
      { Name b ; Value { Term { [ {d a} ]        ; In Domain ; Jacobian Vol ; } } }
      { Name h ; Value { Term { [ nu[{d a}] * {d a} ] ; In Domain ; Jacobian Vol ; } } }

      { Name js; Value { Term { [ js0[] ]      ; In DomainS ; Jacobian Vol ; } } }

      { Name xi ; Value { Term { [ {xi} ]          ; In Domain ; Jacobian Vol ; } } }

      { Name xis ; Value { Term { [ {xis} ] ; In Domain ; Jacobian Vol ; } } }
      { Name dxis ; Value { Term { [ {d xis} ] ; In Domain ; Jacobian Vol ; } } }
      { Name js0_dxis ; Value { Term { [ js0[]-{d xis} ] ; In Domain ; Jacobian Vol ; } } }

      // Magnetic energy and co-energy density, J/m^3
      { Name w_mag ;    Value { Term { [ w_mag[{d a}]    ] ; In Domain ; Jacobian Vol ; } } }
      { Name w_mag_co ; Value { Term { [ w_mag_co[{d a}] ] ; In Domain ; Jacobian Vol ; } } }
      // Magnetic energy and co-energy (J)
      { Name W_mag ; Value {
          Integral { [ SymmetryFactor* w_mag[{d a}] ] ;
            In Domain ; Jacobian Vol ; Integration I_1p ; } } }// Stored in variable $wmag
      { Name W_mag_co ; Value {
          Integral { [ SymmetryFactor* w_mag_co[{d a}] ] ;
            In Domain ; Jacobian Vol ; Integration I_1p ; } } }// Stored in variable $wmag_co

      { Name Flux ; Value {
          Integral { [ SymmetryFactor*vDir[]*NbWires[]/SurfCoil[]*{a} ] ;
            In Inds  ; Jacobian Vol ; Integration II ; } // Stored in variable $flux
        }
      }
      // Inductance computation from flux, energy, co-energy
      { Name L_from_Flux ; Value { Term { Type Global; [ $flux*1e3/II ] ; In DomainDummy ; } } }
      { Name L_from_W_mag; Value { Term { Type Global; [ 2*$wmag*1e3/(II*II) ] ; In DomainDummy ; } } }
      { Name L_from_W_mag_co; Value { Term { Type Global; [ 2*$wmag_co*1e3/(II*II) ] ; In DomainDummy ; } } }


      // Force (N) computation by Maxwell stress tensor (exterior normal automatic from geo, surface orientation matters)
      { Name Fskin  ; Value { Integral { [ SymmetryFactor * Tmax[{d a}] * Normal[] ] ;
            In SkinDomain_Moving ; Jacobian Sur ; Integration I_1p ; } } } // Force_skin_mst stored in variable $fmst_skin
      { Name Fskinx ; Value { Term { Type Global; [ CompX[$fmst_skin] ] ; In DomainDummy ; } } }
      { Name Fskiny ; Value { Term { Type Global; [ CompY[$fmst_skin] ] ; In DomainDummy ; } } }
      { Name Fskinz ; Value { Term { Type Global; [ CompZ[$fmst_skin] ] ; In DomainDummy ; } } }

      // Force computation using Maxwell stress tensor (MST)
      // + Aux quantity for computing correctly the exterior normal
      For i In {1:NumRegionsForceComputation}
        { Name un~{i} ; Value { Local { [ {un~{i}} ] ; In Domain ; Jacobian Vol ; } } }

        { Name Fdensity~{i}  ; Value {
            Local { [ SymmetryFactor*Tmax[{d a}] * {d un~{i}} ] ; In Air ; Jacobian Vol ; } } }
        { Name F~{i}  ; Value { Integral { [ SymmetryFactor*Tmax[{d a}] * {d un~{i}} ] ;
              In Air ; Jacobian Vol ; Integration II ; } } }
        { Name Fx~{i}  ; Value { Integral { [ CompX[SymmetryFactor*Tmax[{d a}] * {d un~{i}}] ] ;
              In Air ; Jacobian Vol ; Integration II ; } } }
        { Name Fy~{i}  ; Value { Integral { [ CompY[SymmetryFactor*Tmax[{d a}] * {d un~{i}}] ] ;
              In Air ; Jacobian Vol ; Integration II ; } } }
        { Name Fz~{i}  ; Value { Integral { [ CompZ[SymmetryFactor*Tmax[{d a}] * {d un~{i}}] ] ;
              In Air ; Jacobian Vol ; Integration II ; } } }
      EndFor

      // checking Normal definion
      { Name N_surf ; Value { Term { [ Normal[] ]    ; In SkinDomain_Moving ; Jacobian Vol ; } } }

    }
  }
}

//-----------------------------------------------------------------------------------------------
 PostOperation Get_LocalFields UsingPost MagSta_a_3D {
   Print[ js, OnElementsOf DomainS, File StrCat[Dir, StrCat["js",ExtGmsh]], LastTimeStepOnly ] ;

   Print[ a, OnElementsOf Domain, File StrCat[Dir, StrCat["a",ExtGmsh]], LastTimeStepOnly ] ;

    If(Flag_DivJ_Zero == DIVJ0_WEAK)
     Print[ xis, OnElementsOf DomainS, File StrCat[Dir, "xis",ExtGmsh ], LastTimeStepOnly ] ;
     Print[ dxis, OnElementsOf DomainS, File StrCat[Dir, "grad_xis",ExtGmsh ], LastTimeStepOnly ] ;
     Print[ js0_dxis, OnElementsOf DomainS, File StrCat[Dir, "js0_corrected",ExtGmsh ], LastTimeStepOnly ] ;
   EndIf

   If(Flag_GaugeType==COULOMB_GAUGE)
     Print[ xi, OnElementsOf Domain, File StrCat[Dir, StrCat["xi",ExtGmsh]], LastTimeStepOnly ] ;
   EndIf
   Print[ b, OnElementsOf Domain, File StrCat[Dir, StrCat["b",ExtGmsh]], LastTimeStepOnly ] ;
 }

 PostOperation Checking_stuff UsingPost MagSta_a_3D {
   Print[ N_surf, OnElementsOf SkinDomain_Moving, File StrCat[Dir,"Nskin",ExtGmsh], LastTimeStepOnly ] ;
 }

 PostOperation Get_GlobalQuantities UsingPost MagSta_a_3D {
   Print[ Flux[DomainS], OnGlobal, Format TimeTable,
     File > StrCat[Dir, StrCat["Flux",ExtGnuplot]], LastTimeStepOnly, StoreInVariable $flux,
     SendToServer StrCat[po,"40Flux [Wb]"],  Color "LightYellow" ];

   Print[ W_mag[Domain], OnGlobal, Format TimeTable,
     File > StrCat[Dir, StrCat["W_mag",ExtGnuplot]], LastTimeStepOnly, StoreInVariable $wmag,
     SendToServer StrCat[po,"41W_mag [J]"],  Color "LightYellow" ];
   Print[ W_mag_co[Domain], OnGlobal, Format TimeTable,
     File > StrCat[Dir, StrCat["W_mag_co",ExtGnuplot]], LastTimeStepOnly, StoreInVariable $wmag_co,
     SendToServer StrCat[po,"42W_mag_co [J]"],  Color "LightYellow" ];

   Print[ L_from_Flux, OnRegion DomainDummy, Format Table, LastTimeStepOnly,
     File StrCat[Dir, StrCat["InductanceF",ExtGnuplot]],
     SendToServer StrCat[poL,"50L from Flux [mH]"], Color "LightYellow" ];
   Print[ L_from_W_mag, OnRegion DomainDummy, Format Table, LastTimeStepOnly,
     File StrCat[Dir, StrCat["InductanceE",ExtGnuplot]],
     SendToServer StrCat[poL,"51L from W_mag [mH]"], Color "LightYellow" ];
   Print[ L_from_W_mag_co, OnRegion DomainDummy, Format Table, LastTimeStepOnly,
     File StrCat[Dir, StrCat["InductanceCoE",ExtGnuplot]],
     SendToServer StrCat[poL,"52L from W_mag_co [mH]"], Color "LightYellow" ];


   Print[ Fskin[SkinDomain_Moving], OnGlobal, Format TimeTable,
     File > StrCat[Dir,StrCat["Fmst_skin",ExtGnuplot]], StoreInVariable $fmst_skin, LastTimeStepOnly ];
   Print[ Fskinx, OnRegion DomainDummy, Format Table, LastTimeStepOnly,
     File  StrCat[Dir,"tmp.dat"], SendToServer StrCat[poF,"60MST (surf)/1Fx [N]"], Color "LightYellow"];
   Print[ Fskiny, OnRegion DomainDummy, Format Table, LastTimeStepOnly,
     File  StrCat[Dir,"tmp.dat"], SendToServer StrCat[poF,"60MST (surf)/2Fy [N]"], Color "NavajoWhite"];
   Print[ Fskinz, OnRegion DomainDummy, Format Table, LastTimeStepOnly,
     File  StrCat[Dir,"tmp.dat"], SendToServer StrCat[poF,"60MST (surf)/3Fz [N]"], Color "LightYellow"];

  // Maxwell stress tensor on a surface, exact computation of normal
  // The explicit definition of a layer around the area of interest is not necessary anymore
  For i In {1:NumRegionsForceComputation}
    // Print[ un~{i}, OnElementsOf Domain, File StrCat[Dir,Sprintf("un_%g.pos",i)] ];
    // Print[ Fdensity~{i}, OnElementsOf Air, File StrCat[Dir,Sprintf("Fd_%g.pos",i)] ];

    Print[ F~{i}[Air],  OnGlobal, Format TimeTable,
      File > StrCat[Dir,Sprintf("Fmst%g.dat",i)], LastTimeStepOnly ];
    Print[ Fx~{i}[Air], OnGlobal, Format Table, LastTimeStepOnly, File StrCat[Dir,"tmp.dat"],
      SendToServer StrCat[poF,"61MST/1Fx [N]"], Color "Ivory"  ];
    Print[ Fy~{i}[Air], OnGlobal, Format Table, File StrCat[Dir,"tmp.dat"],
      SendToServer StrCat[poF,"61MST/2Fy [N]"], Color "NavajoWhite"  ];
    Print[ Fz~{i}[Air], OnGlobal, Format Table, File StrCat[Dir,"tmp.dat"],
      SendToServer StrCat[poF,"61MST/3Fz [N]"], Color "Ivory"  ];
  EndFor
}

DefineConstant[
  R_ = {"Analysis", Name "GetDP/1ResolutionChoices", Visible 1},
  C_ = {"-solve -v 3 -v2 -bin", Name "GetDP/9ComputeCommand", Visible 1},
  P_ = {" ", Name "GetDP/2PostOperationChoices", Visible 1}
];
