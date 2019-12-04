Include "cable_data.geo";

DefineConstant[
  Flag_AnalysisType = {0,
    Choices{
      0="Electric",
      1="Magnetic",
      2="Magneto-thermal (linear)",
      3="Magneto-thermal (nonlinear)"
    },
    Name "{00Parameters/00Type of analysis", Highlight "Blue"}

  Flag_sigma_funcT = (Flag_AnalysisType==3)?1:0,

  nb_iter = 20,
  relaxation_factor = 1,
  stop_criterion = 1e-6,

  r_ = {"Analysis", Name "GetDP/1ResolutionChoices", Visible 1}
  c_ = {"-solve -v2", Name "GetDP/9ComputeCommand", Visible 1},
  p_ = {"", Name "GetDP/2PostOperationChoices", Visible 1, Closed 1}
];

Group {
  /////////////////////////////// Start of lab 3 ////////////////////////////////////
  AirInCable = Region[{AIR_IN}];
  AirAboveSoil = Region[{AIR_OUT}];
  DefectInXLPE = Region[{DEFECT}];
  AirEM = Region[{AIR_IN}];//air in electromagnetic domain
  AirTH = Region[{AIR_OUT}];//air in thermal domain
  Air   = Region[{AIR_OUT,AIR_IN}];//all air

  SemiconductorIn = Region[ {SEMI_IN} ];
  XLPE = Region[{XLPE}];
  SemiconductorOut = Region[ {SEMI_OUT} ];
  APLSheath = Region[{APL}];
  Polyethylene = Region[{POLYETHYLENE_SHEATH,POLYETHYLENE_COVER}];
  SteelA = Region[{STEEL_ARMOUR}];
  SteelP =  Region[{STEEL_PIPE}];
  Steel = Region[{STEEL_ARMOUR,STEEL_PIPE}];

  SoilEM = Region[{SOIL_EM}];//soil in electro-magnetic domain
  SoilTH = Region[{SOIL_TH}];//soil in thermal domain
  Soil  = Region[{SOIL_TH,SOIL_EM}];//all soil

  For k In {1:NbWires}
  Ind~{k} = Region[{(WIRE+k-1)}];//individual conductor
  Inds   += Region[{(WIRE+k-1)}];//all conductors
  EndFor


  Cable = Region[{Inds, SemiconductorIn, XLPE, SemiconductorOut, APLSheath, Polyethylene, Steel, AirInCable}];
  Cable += Region[{DefectInXLPE}];

  Sur_Dirichlet_Ele = Region[{OUTBND_CABLE}];

  // magnetodynamics
  Sur_Dirichlet_Mag = Region[{OUTBND_EM}]; // n.b=0 on this boundary

  DomainCC_Mag  = Region[ {AirEM, Inds} ];
  DomainCC_Mag += Region[ {SemiconductorIn, SemiconductorOut, XLPE, Polyethylene, SoilEM} ];
  DomainC_Mag   = Region[ {Steel, APLSheath} ];

  DomainS0_Mag  = Region[ {} ]; // If imposing source with js0[]
  DomainS_Mag  = Region[ {Inds} ]; // If using Current_2D, it allows accounting for the dependance of sigma with T

  DomainCWithI_Mag  = Region[ {} ];//Leave empty
  Domain_Mag = Region[ {DomainCC_Mag, DomainC_Mag} ];

  // electrodynamics
  Domain_Ele = Region[ {Cable} ]; // Just the cable or the same domain as magnetodynamics

  // Thermal domain
  Vol_Thermal          = Region[{Domain_Mag, AirTH, SoilTH}];
  Vol_QSource_Thermal  = Region[{DomainC_Mag}];
  Vol_QSource0_Thermal = Region[{DomainS0_Mag}];
  Vol_QSourceB_Thermal = Region[{DomainS_Mag}];
  Sur_Convection_Thermal =  Region[{}];
  Sur_Dirichlet_Thermal = Region[{OUTBND_TH}];
  Domain_Thermal = Region[{Vol_Thermal, Sur_Convection_Thermal}];

  DomainDummy = Region[{12345}];
}

Function {
  mu0 = 4.e-7 * Pi;
  eps0 = 8.854187818e-12;

  //permeability
  nu[#{Inds,Soil,Polyethylene,Air}] = 1./mu0;
  nu[Region[{SEMI_IN,SEMI_OUT,XLPE,APL,DEFECT}]]  = 1./mu0;//non ferromagnetic
  nu[Steel]  = 1./(mu0*mur_steel);//ferromagnetic

  // electrical conductivity [S/m]
  sigma[Steel]  = sigma_steel;
  sigma[Polyethylene]  = sigma_polyethylene;
  sigma[Region[{SEMI_IN,SEMI_OUT}] ] = sigma_semiconductor;
  sigma[XLPE] = sigma_xlpe;
  sigma[Soil] = sigma_soil;
  sigma[#{Air,DefectInXLPE}] = 0.;//air

  //resistivity of conductor at temperature $1
  fT_cu[] = (1+alpha_cu*($1-Tref)); // $1 is current temperature in [K], alpha in [1/K]
  fT_al[] = (1+alpha_al*($1-Tref));

  If (!Flag_sigma_funcT)
    sigma[Inds]      = sigma_cu;//non linear
    sigma[APLSheath] = sigma_al;
  Else
    sigma[Inds]      = sigma_cu/fT_cu[$1];//non linear
    sigma[APLSheath] = sigma_al/fT_al[$1];
  EndIf

// relative permittivity
  epsilon[#{Air,Inds,Steel,Soil,APLSheath,DefectInXLPE}] = eps0;
  epsilon[Polyethylene] = eps0*epsr_polyethylene;
  epsilon[Region[{SEMI_IN,SEMI_OUT}]] = eps0*epsr_semiconductor;
  epsilon[Region[{XLPE}]] = eps0*epsr_xlpe;

  Freq = 50;
  Omega = 2*Pi*Freq;

  Pa = 0.; Pb = -120./180.*Pi; Pc = -240./180.*Pi;
  I = 406; // maximum value current in data sheet
  js0[Ind_1] = Vector[0,0,1] * I / SurfaceArea[] * F_Cos_wt_p[]{Omega, Pa};
  js0[Ind_2] = Vector[0,0,1] * I / SurfaceArea[] * F_Cos_wt_p[]{Omega, Pb};
  js0[Ind_3] = Vector[0,0,1] * I / SurfaceArea[] * F_Cos_wt_p[]{Omega, Pc};

  Ns[]= 1;
  Sc[]= SurfaceArea[];

  // second order calculation
  Flag_Degree_a = _deg2_hierarchical ? 2 : 1;
  Flag_Degree_v = _deg2_hierarchical ? 2 : 1;

  // thermal parameters
  Tambient[] = Tamb; // [K]


  // thermal conductivities [W/(m K)]
  k[Steel] = kappa_steel;
  k[APLSheath]  = kappa_al;
  k[Polyethylene]  = kappa_polyethylene;
  k[Region[{SEMI_IN,SEMI_OUT}] ] = kappa_semiconductor;
  k[XLPE] = kappa_xlpe;
  k[Soil] = kappa_soil;
  k[Inds] = kappa_cu;
  k[Air]  = kappa_air;
  // * heat conduction mechanism is the main heat transfer mechanism for an underground cable system
  // * all materials have constant thermal properties, including the thermal resistivity of the soil
  // * radiation and convection are not considered

  // * force convection on ground surface due to wind: h = 7.371 + 6.43*v^0.75
  // Not used here
  h[] = 7.371 + 6.43*v_wind^0.75; // 1, 10 ... Convective coefficient [W/(m^2 K)]
}

Constraint {
  // Electrical constraints
  { Name ElectricScalarPotential;
    Case {
    { Region Ind_1; Value V0; TimeFunction F_Cos_wt_p[]{2*Pi*Freq, Pa}; }
    { Region Ind_2; Value V0; TimeFunction F_Cos_wt_p[]{2*Pi*Freq, Pb}; }
    { Region Ind_3; Value V0; TimeFunction F_Cos_wt_p[]{2*Pi*Freq, Pc}; }

    { Region Sur_Dirichlet_Ele; Value 0;  }
    }
  }
  { Name ZeroElectricScalarPotential; // Only if second order
    Case {
      For k In {1:3}
      { Region Ind~{k}; Value 0; }
    EndFor
    { Region Sur_Dirichlet_Ele; Value 0; }
    }
  }

  // Magnetic constraints
  { Name MagneticVectorPotential_2D;
    Case {
      { Region Sur_Dirichlet_Mag; Value 0; }
    }
  }
  { Name Voltage_2D;
    Case {
    }
  }
  { Name Current_2D;
    Case {
      // constraint used if Inds in DomainS_Mag
      { Region Ind_1; Value I; TimeFunction F_Cos_wt_p[]{2*Pi*Freq, Pa}; }
      { Region Ind_2; Value I; TimeFunction F_Cos_wt_p[]{2*Pi*Freq, Pb}; }
      { Region Ind_3; Value I; TimeFunction F_Cos_wt_p[]{2*Pi*Freq, Pc}; }
    }
  }
/////////////////////////////// End of lab 3 (except Joule losses) ////////////////
/////////////////////////////// Start of lab 4 ////////////////////////////////////

  //////////// Lab 4 thermal ///////////
    // Thermal constraints
   { Name DirichletTemp ;
    Case {
      { Type Assign; Region Sur_Dirichlet_Thermal ; Value Tamb; }
    }
  }

}

Jacobian {
  { Name Vol;
    Case {
      { Region All; Jacobian Vol; }
    }
  }
  { Name Sur;
    Case {
      { Region All; Jacobian Sur; }
    }
  }
}

Integration {
  { Name I1;
    Case {
      { Type Gauss;
        Case {
          { GeoElement Point; NumberOfPoints  1; }
          { GeoElement Line; NumberOfPoints  4; }
          { GeoElement Triangle; NumberOfPoints  4; }
          { GeoElement Quadrangle; NumberOfPoints  4; }

          { GeoElement Triangle2; NumberOfPoints  7; } // Second order, geometrical elements

        }
      }
    }
  }
}

FunctionSpace {

  { Name Hgrad_v_Ele; Type Form0;
    BasisFunction {
      // v = v_n  s_n   ,  for all nodes
      { Name sn; NameOfCoef vn; Function BF_Node;
        Support Domain_Ele; Entity NodesOf[ All ]; }
      If (Flag_Degree_v == 2)
        { Name sn2; NameOfCoef vn2; Function BF_Node_2E;
          Support Domain_Ele; Entity EdgesOf[ All ]; }
      EndIf
    }

    Constraint {
      { NameOfCoef vn; EntityType NodesOf;
        NameOfConstraint ElectricScalarPotential; }
      If (Flag_Degree_v == 2)
        { NameOfCoef vn2;
          EntityType EdgesOf; NameOfConstraint ZeroElectricScalarPotential; }
      EndIf
    }
  }


  { Name Hcurl_a_Mag_2D; Type Form1P;
    BasisFunction {
      { Name se; NameOfCoef ae; Function BF_PerpendicularEdge;
        Support Domain_Mag; Entity NodesOf[ All ]; }
      If (Flag_Degree_a == 2)
        { Name se2; NameOfCoef ae2; Function BF_PerpendicularEdge_2E;
          Support Domain_Mag; Entity EdgesOf[ All ]; }
      EndIf
    }
    Constraint {
      { NameOfCoef ae;
        EntityType NodesOf; NameOfConstraint MagneticVectorPotential_2D; }
      If (Flag_Degree_a == 2)
	{ NameOfCoef ae2; // Only OK if homogeneous BC, otherwise specify zero-BC
          EntityType EdgesOf; NameOfConstraint MagneticVectorPotential_2D; }
      EndIf
    }
  }

  { Name Hregion_i_2D ; Type Vector ;
    BasisFunction {
      { Name sr ; NameOfCoef ir ; Function BF_RegionZ ;
        Support DomainS_Mag ; Entity DomainS_Mag ; }
    }
    GlobalQuantity {
      { Name Is ; Type AliasOf        ; NameOfCoef ir ; }
      { Name Us ; Type AssociatedWith ; NameOfCoef ir ; }
    }
    Constraint {
      { NameOfCoef Us ; EntityType Region ; NameOfConstraint Voltage_2D ; }
      { NameOfCoef Is ; EntityType Region ; NameOfConstraint Current_2D ; }
    }
  }



  // Temperature is discretised with this BF
  { Name Hgrad_Thermal; Type Form0;
    BasisFunction {
      { Name sn; NameOfCoef t; Function BF_Node;
        Support Domain_Thermal; Entity NodesOf[All]; }
    }
    Constraint {

      { NameOfCoef t; EntityType NodesOf ; NameOfConstraint DirichletTemp; }
   }
  }



}

Formulation {

  { Name Electrodynamics_v; Type FemEquation;
    Quantity {
      { Name v; Type Local; NameOfSpace Hgrad_v_Ele; }
    }
    Equation {
      Galerkin { [ sigma[] * Dof{d v} , {d v} ] ;
        In Domain_Ele; Jacobian Vol ; Integration I1 ; }
      Galerkin { DtDof[ epsilon[] * Dof{d v} , {d v} ];
        In Domain_Ele; Jacobian Vol; Integration I1; }
    }
  }

  { Name Darwin_a_2D; Type FemEquation; // Magnetodynamics + displacement current, no coupling
    Quantity {
      { Name a;  Type Local; NameOfSpace Hcurl_a_Mag_2D; }

      { Name ir ; Type Local  ; NameOfSpace Hregion_i_2D ; }
      { Name Us ; Type Global ; NameOfSpace Hregion_i_2D[Us] ; }
      { Name Is ; Type Global ; NameOfSpace Hregion_i_2D[Is] ; }
      { Name T; Type Local ; NameOfSpace Hgrad_Thermal; }
    }
    Equation {
      Galerkin { [ nu[] * Dof{d a} , {d a} ];
        In Domain_Mag; Jacobian Vol; Integration I1; }

      Galerkin { DtDof [ sigma[{T}] * Dof{a} , {a} ];
        In DomainC_Mag; Jacobian Vol; Integration I1; }

      Galerkin { DtDtDof [ epsilon[] * Dof{a} , {a} ]; // Added term => Darwin approximation
        In DomainC_Mag; Jacobian Vol; Integration I1; }

      // Either you impose directly the function js0[]
      Galerkin { [ -js0[] , {a} ];
        In DomainS0_Mag; Jacobian Vol; Integration I1; }

      // or you use the constraints => allows accounting for sigma[{T}]
      Galerkin { [ -Ns[]/Sc[] * Dof{ir}, {a} ] ;
        In DomainS_Mag ; Jacobian Vol ; Integration I1 ; }
      Galerkin { DtDof [ Ns[]/Sc[] * Dof{a}, {ir} ] ;
        In DomainS_Mag ; Jacobian Vol ; Integration I1 ; }
      Galerkin { [ Ns[]/Sc[] / sigma[{T}] * Ns[]/Sc[]* Dof{ir} , {ir} ] ; // resistance term
        In DomainS_Mag ; Jacobian Vol ; Integration I1 ; }
      //GlobalTerm { [ Rdc * Dof{Is} , {Is} ] ; In DomainS ; } // OR this resitance term
      GlobalTerm { [ Dof{Us}, {Is} ] ; In DomainS_Mag ; }

    }
  }



    //Thermal formulation
  { Name ThermalSta ; Type FemEquation;
    Quantity {
      { Name T; Type Local ; NameOfSpace Hgrad_Thermal; }

      // quantities from previous formulations, not unknowns
      { Name a; Type Local ; NameOfSpace Hcurl_a_Mag_2D ; }
      { Name ir ; Type Local  ; NameOfSpace Hregion_i_2D ; }
    }
    Equation {
      Galerkin { [ k[] * Dof{d T} , {d T} ];
	In Vol_Thermal; Integration I1; Jacobian Vol;  }

  //Uncomment for thermal problem

      // Thermal source = Joule losses from EM computation
      // Use <a> when mixing frequency domain and time domain computations
      // {a} is complex, losses are computed taking that into account
    Galerkin { [ -0.5*sigma[{T}]*<a>[SquNorm[Dt[{a}]]], {T}] ;
      In Vol_QSource_Thermal; Integration I1; Jacobian Vol;  }
    /*Galerkin { [ -0.5/sigma[{T}]*<js0[]>[SquNorm[js0[]]], {T} ];
      In Vol_QSource0_Thermal; Integration I1; Jacobian Vol;  }*/
    Galerkin { [ -0.5/sigma[{T}]*<ir>[SquNorm[Ns[]/Sc[]*{ir}]], {T} ];
      In Vol_QSourceB_Thermal; Integration I1; Jacobian Vol;  }

    // Convection boundary condition
    /*Galerkin { [  0] ;
      In Sur_Convection_Thermal; Jacobian Sur ; Integration I1 ; }
    Galerkin { [ 0] ;
      In Sur_Convection_Thermal ; Jacobian Sur ; Integration I1 ; }
    */}
  }
/////////////////////////////////////////////

}

Resolution {

  { Name Analysis;
    System {
      { Name Sys_Ele; NameOfFormulation Electrodynamics_v;
        Type Complex; Frequency Freq; }
      { Name Sys_Mag; NameOfFormulation Darwin_a_2D;
        Type Complex; Frequency Freq; }

        { Name Sys_The; NameOfFormulation ThermalSta; }

    }
    Operation {
      CreateDir["res"];

      If(Flag_AnalysisType == 0) // Electrodynamics
        Generate[Sys_Ele]; Solve[Sys_Ele]; SaveSolution[Sys_Ele];
        PostOperation[Ele_Maps];
        PostOperation[Ele_Cuts];
      EndIf
      If(Flag_AnalysisType == 1) // Magnetodynamics
        InitSolution[Sys_Mag];


        InitSolution[Sys_The]; // Init needed for using the same formulation with sigma dependance or not of {T}

        Generate[Sys_Mag]; Solve[Sys_Mag]; SaveSolution[Sys_Mag];
        PostOperation[Mag_Maps];
        PostOperation[Mag_Global];
      EndIf

      If(Flag_AnalysisType > 1) // Magneto-thermal -- linear or NON linear
      //Uncomment for thermal problem

       InitSolution[Sys_Mag];
       InitSolution[Sys_The];

        If(!Flag_sigma_funcT)
          Generate[Sys_Mag]; Solve[Sys_Mag];
          Generate[Sys_The]; Solve[Sys_The];
        Else
          IterativeLoop[ nb_iter, stop_criterion, relaxation_factor]{
            GenerateJac[Sys_Mag]; SolveJac[Sys_Mag];
            GenerateJac[Sys_The]; SolveJac[Sys_The];
          }
        EndIf
        SaveSolution[Sys_Mag];
        SaveSolution[Sys_The];

        PostOperation[Mag_Maps];
        PostOperation[Mag_Global];
        PostOperation[The_Maps];
      EndIf

    }
  }
}

PostProcessing {

  { Name EleDyn_v; NameOfFormulation Electrodynamics_v;
    Quantity {
      { Name v; Value { Term { [ {v} ]; In Domain_Ele; Jacobian Vol; } } }
      { Name e; Value { Term { [ -{d v} ]; In Domain_Ele; Jacobian Vol; } } }
      { Name em; Value { Term { [ Norm[-{d v}] ]; In Domain_Ele; Jacobian Vol; } } }

      { Name d; Value { Term { [ -epsilon[] * {d v} ]; In Domain_Ele; Jacobian Vol; } } }
      { Name dm; Value { Term { [ Norm[-epsilon[] * {d v}] ]; In Domain_Ele; Jacobian Vol; } } }

      { Name j ; Value { Term { [ -sigma[] * {d v} ] ; In Domain_Ele ; Jacobian Vol; } } }
      { Name jm ; Value { Term { [ Norm[-sigma[] * {d v}] ] ; In Domain_Ele ; Jacobian Vol; } } }

      { Name jtot ; Value {
          Term { [ -sigma[] * {d v} ] ;       In Domain_Ele ; Jacobian Vol; }
          Term { [ -epsilon[] * Dt[{d v}] ] ; In Domain_Ele ; Jacobian Vol; }
        } }

      { Name ElectricEnergy; Value {
          Integral {
            [ 0.5 * epsilon[] * SquNorm[{d v}] ];
            In Domain_Ele; Jacobian Vol; Integration I1;
          }
	}
      }

      { Name V0 ; Value {// For recovering the imposed voltage in post-pro
          Term { Type Global ; [ V0 * F_Cos_wt_p[]{2*Pi*Freq, Pa}] ; In Ind_1 ; }
          Term { Type Global ; [ V0 * F_Cos_wt_p[]{2*Pi*Freq, Pb}] ; In Ind_2 ; }
          Term { Type Global ; [ V0 * F_Cos_wt_p[]{2*Pi*Freq, Pc}] ; In Ind_3 ; }
        } }

      { Name C_from_Energy ; Value { Term { Type Global; [ 2*$We/SquNorm[$voltage]*micro_per_km ] ; In DomainDummy ; } } }
      { Name Ctot ; Value { Term { Type Global; [ 1/(3/$C1) ] ; In DomainDummy ; } } }

      // Stored energy is the same for the 3 phases and the norm[V] too => capacitance per phase is identical
      // Then we have 3 identical capacitances in parallel

    }
  }

  { Name Darwin_a_2D; NameOfFormulation Darwin_a_2D;
    PostQuantity {
      { Name a; Value { Term { [ {a} ]; In Domain_Mag; Jacobian Vol; } } }
      { Name az; Value { Term { [ CompZ[{a}] ]; In Domain_Mag; Jacobian Vol; } } }
      { Name b; Value { Term { [ {d a} ]; In Domain_Mag; Jacobian Vol; } } }
      { Name bm; Value { Term { [ Norm[{d a}] ]; In Domain_Mag; Jacobian Vol; } } }

      //Uncomment for thermal problem

      { Name j; Value { Term { [ -sigma[{T}]*Dt[{a}] ]; In DomainC_Mag; Jacobian Vol; } } }
      { Name jz; Value { Term { [ CompZ[-sigma[{T}]*Dt[{a}]] ]; In DomainC_Mag; Jacobian Vol; } } }
      { Name jm; Value { Term { [ Norm[-sigma[{T}]*Dt[{a}]] ]; In DomainC_Mag; Jacobian Vol; } } }


      { Name d; Value { Term { [ epsilon[] * Dt[Dt[{a}]] ]; In DomainC_Mag; Jacobian Vol; } } }
      { Name dz; Value { Term { [ CompZ[epsilon[] *  Dt[Dt[{a}]] ] ]; In DomainC_Mag; Jacobian Vol; } } }
      { Name dm; Value { Term { [ Norm[epsilon[]  *  Dt[Dt[{a}]] ] ]; In DomainC_Mag; Jacobian Vol; } } }

///////////////////////////////////////////////// Complete fr lab 3///////////////////////////////////////
      { Name JouleLosses; Value {
          Integral { [ 0.5*sigma[{T}]*SquNorm[Dt[{a}]] ]   ; In DomainC_Mag  ; Jacobian Vol ; Integration I1 ; }
          Integral { [ 0.5/sigma[{T}]*SquNorm[js0[]] ]          ; In DomainS0_Mag ; Jacobian Vol ; Integration I1 ; }
          Integral { [ 0.5/sigma[{T}]*SquNorm[Ns[]/Sc[]*{ir}] ] ; In DomainS_Mag  ; Jacobian Vol ; Integration I1 ; }
        }
      }
/////////////////////////////////////////////////////////////////////////////////////////////////////////
      { Name U ; Value {
          Term { [ {Us} ] ; In DomainS_Mag ; }
        }
      }

      { Name I ; Value {
          Term { [ {Is} ] ; In DomainS_Mag ; }
        }
      }

      { Name S ; Value {
          Term { [ {Us}*Conj[{Is}] ] ; In DomainS_Mag ; }
        }
      }

      { Name R ; Value {
          Term { [ -Re[{Us}/{Is}] ] ; In DomainS_Mag ; }
        }
      }

      { Name R_per_km ; Value {
          Term { [ -Re[{Us}/{Is}]*1e3 ] ; In DomainS_Mag ; }
        }
      }

      { Name L ; Value {
          Term { [ -Im[{Us}/{Is}]/(2*Pi*Freq) ] ; In DomainS_Mag ; }
        }
      }

      { Name mL_per_km ; Value {
          Term { [ -1e6*Im[{Us}/{Is}]/(2*Pi*Freq) ] ; In DomainS_Mag ; }
        }
      }

    }
  }

  //Uncomment for thermal problem
    { Name ThermalSta; NameOfFormulation ThermalSta; //NameOfSystem T;
    Quantity {
      { Name T; Value{ Local{ [ {T} ] ; In Vol_Thermal; Jacobian Vol; } } }
      { Name TinC; Value{ Local{ [ {T}-273.15 ] ; In Vol_Thermal; Jacobian Vol; } } }
      { Name q; Value{ Local{ [ -k[]*{d T} ] ; In Vol_Thermal; Jacobian Vol; } } }
    }
  }
}

PostOperation{

  // Electric
  //-------------------------------

  po0 = "{01Capacitance/";

  { Name Ele_Maps; NameOfPostProcessing EleDyn_v;
    Operation {
      Print[ v, OnElementsOf Domain_Ele, File "res/v.pos" ];
      Echo[Str[ // "View[PostProcessing.NbViews-1].ScaleType = 2;", // log
          "View[PostProcessing.NbViews-1].RangeType = 3;" ,// per timestep
          "View[PostProcessing.NbViews-1].ShowTime = 0;",
          "View[PostProcessing.NbViews-1].IntervalsType = 2;"
        ], File "res/v.opt"];
      // Print[ e, OnElementsOf Cable, Name "E [V/m]", File "res/em.pos" ];

      Print[ em, OnElementsOf Cable,
        Name "|E| [V/m]", File "res/em.pos" ];
      Echo[Str[ "View[PostProcessing.NbViews-1].NbIso = 25;",
          "View[PostProcessing.NbViews-1].RangeType = 3;" ,// per timestep
          "View[PostProcessing.NbViews-1].ShowTime = 0;",
          "View[PostProcessing.NbViews-1].IntervalsType = 3;"
        ], File "res/em.opt"];

      Print[ dm, OnElementsOf Cable,
        Name "|D| [A/m²]", File "res/dm.pos" ];
      Echo[Str[ "View[PostProcessing.NbViews-1].NbIso = 25;",
          "View[PostProcessing.NbViews-1].RangeType = 3;" ,// per timestep
          "View[PostProcessing.NbViews-1].ShowTime = 0;",
          "View[PostProcessing.NbViews-1].IntervalsType = 3;"
        ], File "res/dm.opt"];


      Print[ ElectricEnergy[Domain_Ele], OnGlobal, Format Table, StoreInVariable $We,
        SendToServer StrCat[po0,"0Electric energy"], File "res/energy.dat" ];

      Print[ V0, OnRegion Ind_1, Format Table, StoreInVariable $voltage,
        SendToServer StrCat[po0,"0U1"], Units "V", File "res/U.dat" ];
      Print[ C_from_Energy, OnRegion DomainDummy, Format Table, StoreInVariable $C1,
        SendToServer StrCat[po0,"1Cpha"], Units "µF/km", File "res/C.dat" ];

      Print[ Ctot, OnRegion DomainDummy, Format Table,
        SendToServer StrCat[po0,"1Ctot"], Units "µF/km", File >"res/C.dat" ];

    }
  }

  dist_cab = dc + 2*(ti+txlpe+to+tapl)+tps;
  h = dist_cab * Sin[Pi/3]; // height of equilateral triangle
  x0 = 0; y0 = 2*h/3;
  x1 = -dist_cab/2; y1 = -h/3;
  x2 =  dist_cab/2; y2 = -h/3;

  { Name Ele_Cuts; NameOfPostProcessing EleDyn_v;
    Operation {
      Print[ em , OnLine { {x2,y2,0} {x2+dc/2+ti+txlpe+to+tapl,y2,0} } {100},
        Name "|E| [V/m] cut in phase 2", File "res/em_cut.pos"];
      Echo[Str["View[PostProcessing.NbViews-1].Type = 4;",
          "View[PostProcessing.NbViews-1].Axes = 3;",
          "View[PostProcessing.NbViews-1].AutoPosition = 3;", // top right
          "View[PostProcessing.NbViews-1].ShowTime = 0;",
          "View[PostProcessing.NbViews-1].LineWidth = 3;",
          "View[PostProcessing.NbViews-1].NbIso = 5;"
        ],
        File "res/em_cut.opt"];
    }
  }

  // Magnetic
  //-------------------------------

  { Name Mag_Maps; NameOfPostProcessing Darwin_a_2D;
    Operation {
      Print[ bm , OnElementsOf Cable, //Smoothing 1,
        Name "|B| [T]", File "res/bmu.pos" ];
      Echo[Str["View[PostProcessing.NbViews-1].RangeType = 3;" ,// per timestep
          "View[PostProcessing.NbViews-1].ShowTime = 0;",
          "View[PostProcessing.NbViews-1].IntervalsType = 2;"
        ], File "res/maps.opt"];

        //Uncomment for thermal problem

      Print[ jm , OnElementsOf APLSheath ,
        Name "|j| [A/m^2] Al sheath", File "res/jm_al.pos" ];
      Echo[Str["View[PostProcessing.NbViews-1].RangeType = 3;" ,// per timestep
          "View[PostProcessing.NbViews-1].ShowTime = 0;",
          "View[PostProcessing.NbViews-1].IntervalsType = 2;"
        ], File "res/maps.opt"];
      Print[ jm , OnElementsOf SteelA,
        Name "|j| [A/m^2] Steel armour", File "res/jm_steelA.pos" ];
      Echo[Str["View[PostProcessing.NbViews-1].RangeType = 3;" ,// per timestep
          "View[PostProcessing.NbViews-1].ShowTime = 0;",
          "View[PostProcessing.NbViews-1].IntervalsType = 2;"
        ], File "res/maps.opt"];

      Print[ jm , OnElementsOf SteelP,
        Name "|j| [A/m^2] Steel pipe", File "res/jm_steelP.pos" ];
      Echo[Str["View[PostProcessing.NbViews-1].RangeType = 3;" ,// per timestep
          "View[PostProcessing.NbViews-1].ShowTime = 0;",
          "View[PostProcessing.NbViews-1].IntervalsType = 2;"
        ], File "res/maps.opt"];



      Print[ dm , OnElementsOf DomainC_Mag,
        Name "|D| [A/m²]", File "res/dm.pos" ];
      Echo[Str["View[PostProcessing.NbViews-1].RangeType = 3;" ,// per timestep
          "View[PostProcessing.NbViews-1].ShowTime = 0;",
          "View[PostProcessing.NbViews-1].IntervalsType = 2;"
        ], File "res/maps.opt"];
    }
  }

  po = "{01Losses/";
  po2 = "{02PU-parameters/";
  { Name Mag_Global; NameOfPostProcessing Darwin_a_2D;
    Operation {

      Print[ JouleLosses[APLSheath], OnGlobal, Format Table,
        SendToServer StrCat[po,"Aluminum sheath"],
        Units "W/m", File "res/losses_al.dat" ];
      Print[ JouleLosses[SteelA], OnGlobal, Format Table,
        SendToServer StrCat[po,"Steel armour"],
        Units "W/m", File "res/losses_steelA.dat" ];
      Print[ JouleLosses[SteelP], OnGlobal, Format Table,
        SendToServer StrCat[po,"Steel pipe"],
        Units "W/m", File "res/losses_steelP.dat" ];


      Print[ R_per_km, OnRegion Ind_1, Format Table,
        SendToServer StrCat[po2,"0R"],
        Units "Ω/km", File "res/Rinds.dat" ];
      Print[ mL_per_km, OnRegion Ind_1, Format Table,
        SendToServer StrCat[po2,"1L"],
        Units "m/km", File "res/Linds.dat" ];

    }
  }

  //Uncomment for thermal problem

  // Thermal
  // -------------------------------
  { Name The_Maps; NameOfPostProcessing ThermalSta;
    Operation {
      Print[ TinC , OnElementsOf Region[{Vol_Thermal,-Cable}], Name "T [°C] around cable", File "res/T.pos" ];
      Echo[Str["View[PostProcessing.NbViews-1].RangeType = 3;" ,// per timestep
          "View[PostProcessing.NbViews-1].ShowTime = 0;",
          "View[PostProcessing.NbViews-1].IntervalsType = 3;",
          "View[PostProcessing.NbViews-1].NbIso = 10;",
          "View[PostProcessing.NbViews-1].ColormapNumber = 6;",
          // "View[PostProcessing.NbViews-1].ColormapSwap = 1;",
          "View[PostProcessing.NbViews-1].Axes = 2;" // Axes (0: none, 1: simple axes, 2: box, 3: full grid, 4: open grid, 5: ruler)
        ], File "res/maps.opt"];
      Print[ TinC , OnElementsOf Cable, Name "T [°C] cable", File "res/Tcable.pos" ];
      Echo[Str["View[PostProcessing.NbViews-1].RangeType = 3;" ,// per timestep
          "View[PostProcessing.NbViews-1].ShowTime = 0;",
               "View[PostProcessing.NbViews-1].IntervalsType = 2;",
               "View[PostProcessing.NbViews-1].ColormapNumber = 6;"
        ], File "res/maps.opt"];
    }
  }

}
