Include "inductor3d_data.geo";

SetFactory("OpenCASCADE");

Mesh.Optimize = 1;

DefineConstant[
  md = { 1.,  Name StrCat[ppm,"0Mesh density"], Highlight Str[colorpp], Closed close_menu},
  nn_wcore   = { Ceil[md*4],  Name StrCat[ppm,"0core width"], ReadOnly 1, Highlight Str[colorro], Closed close_menu},
  nn_airgap  = { Ceil[md*1], Name StrCat[ppm,"1air gap width"], ReadOnly 1, Highlight Str[colorro]},
  nn_ri = { Ceil[md*12], Name StrCat[ppm,"2"], Label "1/4 shell in", ReadOnly 1, Visible (Flag_Infinity==1), Highlight Str[colorro]},
  nn_ro = { Ceil[md*12], Name StrCat[ppm,"3"], Label "1/4 shell out", ReadOnly 1, Highlight Str[colorro]}
];

// characteristic lengths
lc0  = wcoil/nn_wcore;
lc1  = ag/nn_airgap;
lc2  = 4*lc1;
lcri = Pi*Rint/2/nn_ri;

// -----------------------------------------------------------------------------
// E-core
xe = -wcoreI/2*0;          dxe = wcoreI/2;
ye = -htot/2+wcoreE+ag ; dye = hcoreE;
ze = -Lz/2;              dze = Lz/2;
vE()+=newv; Box(newv) = {xe, ye,  ze, dxe, dye, dze};

// I-core
xi = -wcoreI/2*0; dxi = wcoreI/2;
yi = -htot/2 ;  dyi = wcoreE;
zi = -Lz/2; dzi = Lz/2;
vCoreI()+=newv; Box(newv) = {xi, yi, zi, dxi, dyi, dzi};

// coil
xc = -wcoreI/2*0+wcoreE; dxc = wcoil;
yc = -htot/2+wcoreE+ag ; dyc = hcoil;
zc = -Lz/2;              dzc = Lz/2;
vBc()+=newv; Box(newv) = {xc, yc,  zc, dxc, dyc, dzc}; // part of coil intersecting the E-Core
vBc()+=newv; Box(newv) = {0., yc,  zc, wcoreE, dyc, -wcoil}; // part of coil outside the E-Core
vCc()+=newv; Cylinder(newv) = {0, -htot/2+wcoreE+ag, 0, 0,  hcoil, 0, wcoil, Pi/2}; // bend of coil
Rotate {{0, 1, 0}, {0, 0, 0}, Pi/2} { Volume{vCc(0)}; }
Translate {wcoreE, 0,  -Lz/2} { Volume{vCc(0)};} // also touching the E-Core

vCoil()+=newv; BooleanUnion(newv) = { Volume{vBc()}; Delete; }{ Volume{vCc()}; Delete; }; // Complete coil

vCoreE() += newv; BooleanDifference(newv) = { Volume{vE(0)}; Delete; }{ Volume{vCoil()}; };

// Air around
vA()+=newv; Sphere(newv) = {0,0,0, Rint, -Pi/2, Pi/2, Pi/2};
vA()+=newv; Sphere(newv) = {0,0,0, Rext, -Pi/2, Pi/2, Pi/2};
Rotate {{1, 0, 0}, {0, 0, 0}, -Pi/2} { Volume{vA()}; }

vAir()+=newv; BooleanDifference(newv) = { Volume{vA(1)}; Delete; }{ Volume{vA(0)}; };
vAir()+=newv; BooleanDifference(newv) = { Volume{vA(0)}; Delete; }{ Volume{vCoreE(),vCoreI(),vCoil()}; };

BooleanFragments{ Volume{vAir(), vCoreE(), vCoreI(), vCoil()}; Delete; }{} // This needs to be done at the end
// Coherence; // This does exactly the same as BooleanFragments when what we want is to paste parts of the geometry

// Adapting mesh size...
Characteristic Length { PointsOf{ Volume{vAir(0)}; } }= lcri;
Characteristic Length { PointsOf{ Volume{vCoreE(), vCoreI()}; } }= lc0;
Characteristic Length { PointsOf{ Volume{vCoil()}; } }= lc2; // Basic lc, we refine after

// and around the airgap
ptos_ag() = Point In BoundingBox  {xi-ag, -htot/2+wcoreE-ag, zi-ag, dxi, 3*ag, dzi};
Characteristic Length { ptos_ag() }= lc1; // points around airgap


// Boundary conditions
tol = 2*ag;
cut_xy() = Surface In BoundingBox {-Rext-tol,-Rext-tol,-tol, 2*(Rext+tol), 2*(Rext+tol), 2*tol}; // 1/2, 1/4
cut_yz() = Surface In BoundingBox {-tol,-Rext-tol,-Rext-tol, 2*tol, 2*(Rext+tol), 2*(Rext+tol)}; // 1/4

all_vol() = Volume{:};
bndAir() = CombinedBoundary{Volume{all_vol()};};
bndAir() -= {cut_xy(),cut_yz()};


//=================================================
// Some colors... for aesthetics :-)
//=================================================
Recursive Color SkyBlue {Volume{vAir()};}
Recursive Color SteelBlue {Volume{vCoreE(),vCoreI()};}
Recursive Color Red {Volume{vCoil()};}

//=================================================
// Physical regions for FE analysis with GetDP
//=================================================

Physical Volume("E core", ECORE) = vCoreE();
Physical Volume("I core", ICORE) = vCoreI();

bnd_vCoreI() = Boundary{Volume{vCoreI()};};
bnd_vCoreI() -= {cut_xy(),cut_yz()};
Physical Surface("Boundary I core", SKINICORE) = bnd_vCoreI();

bnd_vCoreE() = CombinedBoundary{Volume{vCoreE()};};
bnd_vCoreE() -= {cut_xy(),cut_yz()};
Physical Surface("Boundary E core", SKINECORE) = bnd_vCoreE();

Physical Volume(COIL) = vCoil();
bnd_vCoil() = CombinedBoundary{Volume{vCoil()};};
bnd_vCoil() -= {cut_xy(),cut_yz()};
Physical Surface("Boundary Coil", SKINCOIL) = bnd_vCoil();

If(Flag_Infinity==0)
Physical Volume("Air", AIR) = {vAir()};
EndIf
If(Flag_Infinity==1)
  Physical Volume("Air", AIR) = vAir(1);
  Physical Volume("Air infinity shell", AIRINF) = vAir(0);
EndIf
Physical Surface("Outer boundary", SURF_AIROUT) = bndAir();

Physical Surface("Symmetry cut xy", CUT_XY) = {cut_xy()};
Physical Surface("Symmetry cut yz", CUT_YZ) = {cut_yz()}; // BC if symmetry
