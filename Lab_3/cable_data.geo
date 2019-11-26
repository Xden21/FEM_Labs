//=================================================
// Geometrical data
//=================================================

mm = 1e-3;

dc = 18.4*mm; // Diameter of conductor
ti = 1.4 *mm;  // Thickness of inner semi-conductor
txlpe = 11*mm; // Thickness of XLPE insulation

to   = 1.4*mm; // Thichness of outher semi-conductor
tapl = 2*mm; // Thickness of APL sheath
tps  = 1*mm;  // Thickness of Polyethylene sheath
tswa = 2*2*mm; // Thickness of Steel wire armour (+++ I think there is an error in the Table I, geometrical data)
tsp  = 4*mm;  // Thickness of Steel pipe
tpc  = 1*mm;  // Thickness of Polyethylene covering

depth_cable=1.2; // [m] laying depth of the cable

dtot = 135*mm; // Overall diameter dc
dinf = 5*dtot; // Electromagnetic analysis

// dinf_th = 8*dinf; // thermal analysis
// dinf_th_air = 10*dtot;

dinf_th = 3; // thermal analysis
dinf_th_air = 3 -depth_cable;



DefineConstant[
  _use_2ndorder_geo = {0, Choices{0,1},
    Name "{00Parameters/000Second order geometrical elements? (experimental)", Highlight "LightYellow"}
  _deg2_hierarchical = {(_use_2ndorder_geo?0:1), Choices{0,1}, Visible 1, ReadOnly 1,
    Name "{00Parameters/001Second order hierarchical elements?"}

  Flag_defect_in_XLPE = {0,  Choices{0,1},
    Name "{00Parameters/10Defect in XLPE (phase 2)", Highlight "Red"}
  dd = {2, Choices{2,5},Name "Parameters/11Distance from conductor defect", Visible Flag_defect_in_XLPE}
];
dd = dd*mm; // position of defect
ddefect = 2*mm; // => defect diameter in XLPE layer, phase 2 (from figures)

NbWires=3;

//=================================================
// material properties
//=================================================
// relative permittivity
epsr_polyethylene = 2.25;
epsr_semiconductor = 2.25;
epsr_xlpe = 2.5;
// = 1 for Cu, steel, Al, soil (dry)

// relative permeability
mur_steel = 4;
// 1 for Cu, Al, polyethylene, semiconductor, XLPE, soil(dry)

// electrical conductivity [S/m]
sigma_cu    = 5.99e7;
sigma_steel = 4.7e6;
sigma_al    = 3.77e7;
sigma_polyethylene =  1.0e-18;
sigma_semiconductor = 2;
sigma_xlpe = 1.0e-18;
sigma_soil = 28;

// thermal conductivity [W/(m K)]
kappa_cu   = 400 ;
kappa_steel= 50.2;
kappa_al   = 237;
kappa_polyethylene   = 0.46;
kappa_semiconductor = 10;
kappa_xlpe = 0.46;
kappa_soil = 0.4;

kappa_air = 0.025499; // typical air at Tamb = 15 and 1 atm (from IPSE course)

//=================================================
// Parameters for simulations
//=================================================

// electrical analysis
Vrms = 132e3; // RMS value in the line voltage [V]
V0 = Vrms/Sqrt[3]; // peak value
// add a defect in the XLPE insulation => d=5mm or d=2mm

// magnetic analysis
// coupling with thermal problem
// conductivity depends on temperature sigma(T)
Tref = 20+273.15; // [K] ==20 C

// alpha is equal to 0.00386K-1 for the copper and 0.00390K-1 for the aluminum
// sigma(T) =  sigma0/(1+alpha*(T-Tref));
alpha_cu = 0.00386; //[1/K]
alpha_al = 0.00390; //[1/K]
// conductivity of steel does not vary with temperature


// The ambient temperature considered in the thermal analysis is T=15°C.
// The iterative non-linear solver found an accurate solution in five iterations and about 25s.
// The maximum temperature of the cable is 77.35 °C, that is the operating temperature of real cables
// (about 70 to 90°C).
// When the cable carries the maximum current value the soil surface temperature raise up from 15°C to about 35°C.
Tamb = 15+273.15; // [K] ==15 C, ambient temperature for the thermal analysis

v_wind = 1.2; // [m/s] wind speed for force convection at interface air/ground

eps0 = 8.854187818e-12;
R1 = dc/2;
R2 = dc/2+2*ti+txlpe;
Cana = 2*Pi*eps0*epsr_xlpe/Log[R2/R1];
micro_per_km = 1e3/1e-6;

Printf("===> Cana %g µF/km",Cana*micro_per_km);


//=================================================
// Physical numbers
//=================================================
AIR_IN = 900; AIR_OUT=901;

WIRE  = 1000;

SEMI_IN = 2000;
XLPE = 3000; DEFECT=3001;
SEMI_OUT =4000;
APL = 5000;
POLYETHYLENE_SHEATH = 6000;

STEEL_ARMOUR = 7000;
STEEL_PIPE = 8000;

POLYETHYLENE_COVER = 9000;
SOIL_EM = 10000;
SOIL_TH = 11000;

OUTBND_EM=1111;
OUTBND_TH=2222;

INTERFACE_AIR_SOIL=3333;

OUTBND_CABLE=1110;
