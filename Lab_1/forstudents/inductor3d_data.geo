// Geometrical data for inductor model

mm = 1e-3; // Unit

pp  = '{1Geometrical data/10dimensions/0';
pp2 = '{1Geometrical data/10dimensions/1Shell radius/';
ppm = '{1Geometrical data/11mesh control/0';

psim= '{2Simulation param./';

DefineConstant[
  Flag_Symmetry = {2, Choices{0="Full",1="Half",2="One fourth"},
    Name '{1Geometrical data/00Symmetry type', Highlight 'Blue', ReadOnly 1}
  Flag_Infinity = {1, Choices{0,1},
    Name '{1Geometrical data/01Use shell transformation to infinity', Highlight 'White', Visible 1}
  Flag_OpenCore = 1
];

SymmetryFactor = 4 ;

close_menu = 1;
colorro  = "LightGrey";
colorpp = "Ivory";

DefineConstant[
  wI = {25.4,  Name StrCat[pp,'01total E-core width'], Units 'mm' ,Highlight Str[colorpp],Closed close_menu}
  wcoreE = {wI/8, ReadOnly 1, Name StrCat[pp,'02E-core width of side legs'], Units 'mm', Highlight Str[colorro]}
  hcoil  = {wI/2, ReadOnly 1, Name StrCat[pp,'04Coil height'], Units 'mm', Highlight Str[colorpp]},
  wcoil  = {wI/4, ReadOnly 1, Name StrCat[pp,'03Coil width'], Units 'mm', Highlight Str[colorro]},
  hcoreE = {hcoil+wcoreE, ReadOnly 1, Name StrCat[pp,'02E-core height of legs'], Units 'mm', Highlight Str[colorro]},
  ag     = { wcoreE/10, Min 0., Max wcoreE/5, Step (wcoreE/5)/10, Units 'mm', ReadOnlyRange 1,
    Name StrCat[pp,'5Air gap width'], Highlight Str[colorpp]},
  Lz = {18*0.36, Name StrCat[pp,'09Core length along z-axis'], Units 'mm', Highlight Str[colorpp]}
];

// laminated core
// 18 lams x 0.36 mm = 6.48 mm
// 14 lams x 0.47 mm = 6.58 mm

// Converting to meters
wI = wI*mm;
wcoreE = wcoreE*mm;
hcoil  = hcoil*mm;
wcoil  = wcoil*mm;
hcoreE = hcoreE*mm;
ag = ag*mm;
Lz = Lz*mm;


// rest of EI-Core dimensions
wcoreE_centralleg = 2*wcoreE;

wcoreI = wI;
hcoreI = wcoreE;

htot = hcoreE + ag + wcoreE ; // Total height of EI-core, including gap



alw = ag/3; // width of air layer around I-core for force computation

// radious for surrounding air with transformation to infinity

If(Flag_Infinity==1)
  label_Rext = "1Outer [m]";
EndIf
If(Flag_Infinity==0)
  label_Rext = "1[m]";
EndIf

DefineConstant[
  Rint = {30*mm, Min 0.15, Max 0.9, Step 0.1,
    Name StrCat[pp2,'0Inner [m]'], Visible (Flag_Infinity==1), Highlight Str[colorpp] },
  Rext = {40*mm, Min Rint, Max 1, Step 0.1,
    Name StrCat[pp2, label_Rext], Visible 1, Highlight Str[colorpp] }
];

Val_Rint = Rint;
Val_Rext = Rext;

IA = 5;
Nw = 60;

sigma_al = 3.72e7 ; // conductivity of aluminum [S/m]
sigma_cu = 5.77e7  ; // conductivity of copper [S/m]

// ----------------------------------------------------
// Numbers for physical regions in .geo and .pro files
// ----------------------------------------------------

ECORE = 1000;
ICORE = 1100;

SKINECORE = 1111;
SKINICORE = 1112;

COIL  = 2000;
LEG_INCOIL = 2100;
SKINCOIL = 2222;
SKINCOIL_ = 2223;

SURF_ELEC0 = 2333;
SURF_ELEC1 = 2444;
CUTCOIL = 2555;

AIR    = 3000;
AIRINF = 3100;
AIRGAP = 3200;

AIRLAYER = 3300;

AIRLAYER_U = 3301;
AIRLAYER_D = 3302;
AIRLAYER_R = 3303;
AIRLAYER_B = 3304;

// Lines and surfaces for boundary conditions
SURF_AIROUT = 3333;

AXIS_Y = 10000;
CUT_YZ = 11000;
CUT_XY = 12000;
