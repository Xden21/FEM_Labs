Include "inductor2d_data.geo";

// Gmsh project created on Wed Nov 06 17:44:14 2019
SetFactory("OpenCASCADE");

// characteristic lengths (length of side of 1 mesh triangle)
lc0  = wC/nn_wcore; // Used for core
lc1  = ag/nn_airgap; // Used for gap
lc2  = 4*lc1;  // Used for turn at top
lcri = Pi*Rint/2/nn_ri; // Used for border cirkles

// Define center
center = newp; Point(newp) = {0,0,0,lc0};

// E-Core
// Define bottom line with 4 points (in array for convenience) use newp for next available point number
// CLose to air gap so small mesh
epbtm[] += newp; Point(newp) = {0, htot/2 - hE, 0, lc1};
epbtm[] += newp; Point(newp) = {wE, htot/2 - hE, 0, lc1};
epbtm[] += newp; Point(newp) = {wE + wC, htot/2 - hE, 0, lc1};
epbtm[] += newp; Point(newp) = {2*wE + wC, htot/2 - hE, 0, lc1};

// Define Top line E-Core with 2 points
eptop[] += newp; Point(newp) = {0, htot/2, 0, lc0};
eptop[] += newp; Point(newp) = {2*wE + wC, htot/2, 0, lc0};

// Define Top coil line with 2 points
epcoil[] += newp; Point(newp) = {wE, htot/2-hE+hC, 0, lc2};
epcoil[] += newp; Point(newp) = {wE + wC, htot/2-hE+hC, 0, lc2};

// Draw E-core lines

// Horizontal
elbtm[] += newl; Line(newl) = {epbtm[0], epbtm[1]};
elbtm[] += newl; Line(newl) = {epbtm[1], epbtm[2]};
elbtm[] += newl; Line(newl) = {epbtm[2], epbtm[3]};

eltop = newl; Line(newl) = {eptop[0], eptop[1]};
elcoil = newl; Line(newl) = {epcoil[0], epcoil[1]};

// Vertical
elvert[] += newl; Line(newl) = {epbtm[0], eptop[0]};
elvert[] += newl; Line(newl) = {epbtm[1], epcoil[0]};
elvert[] += newl; Line(newl) = {epbtm[2], epcoil[1]};
elvert[] += newl; Line(newl) = {epbtm[3], eptop[1]};

// I-Core

// Define corner points
ipbtm[] += newp; Point(newp) = {0 , -htot/2, 0, lc0};
ipbtm[] += newp; Point(newp) = {wI/2 , -htot/2, 0, lc0};

iptop[] += newp; Point(newp) = {0, -htot/2 + hI, 0, lc1};
iptop[] += newp; Point(newp) = {wI/2, -htot/2 + hI, 0, lc1};

// Define border lines
// Horizontal
 ilhor[] += newl; Line(newl) = {iptop[0], iptop[1]};
 ilhor[] += newl; Line(newl) = {ipbtm[0], ipbtm[1]};

 // Vertical
 ilvert[] += newl; Line(newl) = {iptop[0], ipbtm[0]};
 ilvert[] += newl; Line(newl) = {iptop[1], ipbtm[1]};

 // Connect points of airgap
 airvert[] += newl; Line(newl) = {epbtm[0], iptop[0]};
 airvert[] += newl; Line(newl) = {epbtm[3], iptop[1]};

// Inner Circle
// Define 3 points
rinpts[] += newp; Point(newp) = {0, Rint, 0, lcri};
rinpts[] += newp; Point(newp) = {0, -Rint, 0, lcri};
rinpts[] += newp; Point(newp) = {Rint, 0, 0, lcri};

// Define circles
lnrin[] += newl; Circle(newl) = {rinpts[0], center, rinpts[2]};
lnrin[] += newl; Circle(newl) = {rinpts[2], center, rinpts[1]};

// close lines
lnrinlines[] += newl; Line(newl) = {rinpts[0], eptop[0]};
lnrinlines[] += newl; Line(newl) = {rinpts[1], ipbtm[0]};

// Outer Circle
// Define 3 points
rextpts[] += newp; Point(newp) = {0, Rext, 0, lcri};
rextpts[] += newp; Point(newp) = {0, -Rext, 0, lcri};
rextpts[] += newp; Point(newp) = {Rext, 0, 0, lcri};

// Define circles
lnrext[] += newl; Circle(newl) = {rextpts[0], center, rextpts[2]};
lnrext[] += newl; Circle(newl) = {rextpts[2], center, rextpts[1]};

// close lines
lnrextlines[] += newl; Line(newl) = {rextpts[0], rinpts[0]};
lnrextlines[] += newl; Line(newl) = {rextpts[1], rinpts[1]};

// E Core
Curve Loop(1) = {6, 4, -9, -3, 8, -5, -7, -1};
//+
Surf_E_Core = news;  Plane Surface(news) = {1};
//+
// Coil
Curve Loop(2) = {2, 8, -5, -7};
//+
Surf_Coil = news; Plane Surface(news) = {2};
//+
// I Core
Curve Loop(3) = {10, 13, -11, -12};
//+
Surf_I_Core = news; Plane Surface(news) = {3};
//+
// Air Gap
Curve Loop(4) = {10, -15, -3, -2, -1, 14};
//+
Surf_Air_Gap = news; Plane Surface(news) = {4};

// Air around inductor inside Rin
Curve Loop(5) = {18, 4, -9, 15, 13, -11, -19, -17, -16};
//+
Surf_Air_Rin = news; Plane Surface(news) = {5};

// Air Between Rin and Rext
Curve Loop(6) = {22, 16, 17, -23, -21, -20};
//+
Surf_Air_Rext = news; Plane Surface(news) = {6};

// Physical surfaces
//+
Physical Surface("E-core", ECORE) = Surf_E_Core[];
//+
Physical Surface("I-core", ICORE) = Surf_I_Core[];
//+
Physical Surface("Coil", COIL) = Surf_Coil[];

Physical Surface("Air", AIR) = {Surf_Air_Rin[],Surf_Air_Gap[]};

Physical Surface("Air Inf", AIRINF) = Surf_Air_Rext[];

Physical Surface("Air Gap", AIRGAP) = Surf_Air_Gap[];
//+
Physical Curve("E-core Skin", SKINECORE) = {elbtm[0], elbtm[2], eltop[], elvert[1], elvert[2], elvert[3], elcoil[]};
//+
Physical Curve("I-core Skin", SKINICORE) = {ilhor[], ilvert[1]};

Physical Curve("Y-Axis", AXIS_Y) = {lnrinlines[], lnrextlines[], elvert[0], ilvert[0], airvert[0]};

Physical Curve("Coil Skin", SKINCOIL) = {elbtm[1], elvert[1], elvert[2], elcoil[]};

Physical Curve("Skin Core and Coil", SKINCORE_COIL) = {elvert[3], elbtm[], eltop[]};

Physical Curve("Outer Edge", SURF_AIROUT) = lnrext[];
