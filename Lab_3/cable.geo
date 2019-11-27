Include "cable_data.geo";

Mesh.Algorithm = 6;
// 2D mesh algorithm (1: MeshAdapt, 2: Automatic, 5: Delaunay, 6: Frontal-Delaunay, 7: BAMG, 8: Frontal-Delaunay for Quads, 9: Packing of Parallelograms)

Mesh.ElementOrder = (!_use_2ndorder_geo)?1:2;

SetFactory("OpenCASCADE");

dist_cab = dc + 2*(ti+txlpe+to+tapl)+tps;
h = dist_cab * Sin(Pi/3); // height of equilateral triangle

x0 = 0; y0 = 2*h/3;
x1 = -dist_cab/2; y1 = -h/3;
x2 =  dist_cab/2; y2 = -h/3;

sur_wire()={};
sur_wire(0) = news; Disk(news) = {x0, y0, 0., dc/2};
sur_wire(1) = news; Disk(news) = {x1, y1, 0., dc/2};
sur_wire(2) = news; Disk(news) = {x2, y2, 0., dc/2};

sur_semi_in()={};
sur_semi_in(0) = news; Disk(news) = {x0, y0, 0., dc/2+ti};
sur_semi_in(1) = news; Disk(news) = {x1, y1, 0., dc/2+ti};
sur_semi_in(2) = news; Disk(news) = {x2, y2, 0., dc/2+ti};

sur_xlpe()={};
sur_xlpe(0) = news; Disk(news) = {x0, y0, 0., dc/2+ti+txlpe};
sur_xlpe(1) = news; Disk(news) = {x1, y1, 0., dc/2+ti+txlpe};
sur_xlpe(2) = news; Disk(news) = {x2, y2, 0., dc/2+ti+txlpe};

sur_semi_out()={};
sur_semi_out(0) = news; Disk(news) = {x0, y0, 0., dc/2+ti+txlpe+to};
sur_semi_out(1) = news; Disk(news) = {x1, y1, 0., dc/2+ti+txlpe+to};
sur_semi_out(2) = news; Disk(news) = {x2, y2, 0., dc/2+ti+txlpe+to};

sur_al()={};
sur_al(0) = news; Disk(news) = {x0, y0, 0., dc/2+ti+txlpe+to+tapl};
sur_al(1) = news; Disk(news) = {x1, y1, 0., dc/2+ti+txlpe+to+tapl};
sur_al(2) = news; Disk(news) = {x2, y2, 0., dc/2+ti+txlpe+to+tapl};


// rounded triangles
R0 = dc/2+ti+txlpe+to+tapl+0*tps;
R1 = R0 + tps;
R2 = R1 + tswa;

arc_angle = 2*Pi/3;

arc_aux0() += newl; Circle(newl) = {x0,y0, 0, R0, Pi/6,     Pi/6+arc_angle};
arc_aux0() += newl; Circle(newl) = {x1,y1, 0, R0, 5*Pi/6, 5*Pi/6+arc_angle};
arc_aux0() += newl; Circle(newl) = {x2,y2, 0, R0, 9*Pi/6, 9*Pi/6+arc_angle};
lin_aux() += newl; Line(newl) = {21, 18};
lin_aux() += newl; Line(newl) = {17, 20};
lin_aux() += newl; Line(newl) = {19, 16};
Curve Loop(newll) = {21, 16, 20, 18, 19, 17};
sur_tri()+=news; Plane Surface(news) = {newll-1};

arc_aux1() += newl; Circle(newl) = {x0,y0, 0, R1, Pi/6,     Pi/6+arc_angle};
arc_aux1() += newl; Circle(newl) = {x1,y1, 0, R1, 5*Pi/6, 5*Pi/6+arc_angle};
arc_aux1() += newl; Circle(newl) = {x2,y2, 0, R1, 9*Pi/6, 9*Pi/6+arc_angle};
lin_aux() += newl; Line(newl) = {23, 26};
lin_aux() += newl; Line(newl) = {27, 24};
lin_aux() += newl; Line(newl) = {25, 22};
Curve Loop(newll) = {27, 26, 28, 25, 29, 24};
sur_tri()+=news; Plane Surface(news) = {newll-1};

arc_aux2() += newl; Circle(newl) = {x0,y0, 0, R2, Pi/6,     Pi/6+arc_angle};
arc_aux2() += newl; Circle(newl) = {x1,y1, 0, R2, 5*Pi/6, 5*Pi/6+arc_angle};
arc_aux2() += newl; Circle(newl) = {x2,y2, 0, R2, 9*Pi/6, 9*Pi/6+arc_angle};
lin_aux() += newl; Line(newl) = {29, 32};
lin_aux() += newl; Line(newl) = {33, 30};
lin_aux() += newl; Line(newl) = {31, 28};
Curve Loop(newll) = {32, 35, 34, 36, 33, 37};
sur_tri()+=news; Plane Surface(news) = {newll-1};

sur_covering()={};
sur_covering(0)=news; Disk(sur_covering(0)) = {0., 0., 0., dtot/2};
sur_covering(1)=news; Disk(sur_covering(1)) = {0., 0., 0., dtot/2-tpc};
sur_covering(2)=news; Disk(sur_covering(2)) = {0., 0., 0., dtot/2-tpc-tsp};

BooleanFragments{
  Surface{
    sur_wire(), sur_semi_in(), sur_xlpe(), sur_semi_out(), sur_al(),
    sur_tri(),
    sur_covering()}; Delete; }{ }

sur_wire() = {1,2,3};
sur_semi_in() = {4,5,6};
sur_xlpe() = {7,8,9};
sur_semi_out() = {10,11,12};
sur_al() = {13,14,15};
sur_ps = 17;
sur_steel_armour=18;
sur_air() = {20,16};
sur_steel_pipe=21;
sur_pc = 19;


// Adding further air detail in cable
sur_ps_aux()={};
sur_ps_aux(0) = news; Disk(news) = {x0, y0, 0., R0+tps};
sur_ps_aux(1) = news; Disk(news) = {x1, y1, 0., R0+tps};
sur_ps_aux(2) = news; Disk(news) = {x2, y2, 0., R0+tps};

sur_wire_cover0() = {sur_wire(0),sur_semi_in(0),sur_xlpe(0),sur_semi_out(0),sur_al(0)};
sur_wire_cover1() = {sur_wire(1),sur_semi_in(1),sur_xlpe(1),sur_semi_out(1),sur_al(1)};
sur_wire_cover2() = {sur_wire(2),sur_semi_in(2),sur_xlpe(2),sur_semi_out(2),sur_al(2)};
sur_ps2(0)=news; BooleanDifference(news) = { Surface{sur_ps_aux(0)}; Delete; }{ Surface{sur_wire_cover0()}; };
sur_ps2(1)=news; BooleanDifference(news) = { Surface{sur_ps_aux(1)}; Delete; }{ Surface{sur_wire_cover1()}; };
sur_ps2(2)=news; BooleanDifference(news) = { Surface{sur_ps_aux(2)}; Delete; }{ Surface{sur_wire_cover2()}; };

sur_air_in()   = BooleanDifference{ Surface{sur_air(1)}; Delete;}{ Surface{sur_ps2()}; };
sur_air()-=sur_air(1);
sur_ps_exact() = BooleanFragments{Surface{sur_ps2(),sur_ps}; Delete; }{ };

/////////////////////////////// to complete in Lab 3 /////////////////////////////////////
//Introduce a defect in phase 2 of the XLPE insulation
If(Flag_defect_in_XLPE)
//Hint: use BooleanDifference
  sur_defect = news;
  Disk(news) = *******; // ti insulation thickness removed
  sur_xlpe(2)= *******;
EndIf
/////////////////////////////////////////////////////////////////////////////////////////

//----------------------------------
// Around the cable
//----------------------------------

all_sur_cable() = Surface{:};

// electromagnetic analysis domain
sur_EMdom = news; Disk(news) = {0, 0, 0., dinf/2};
BooleanDifference(news) = { Surface{sur_EMdom}; Delete; }{ Surface{all_sur_cable()}; };
sur_EMdom=news-1;



// thermal analysis domain
sur_soil   = news; Rectangle(news) = {-dinf_th, -dinf_th, 0, 2*dinf_th, dinf_th+depth_cable, 0};
BooleanDifference(news) = { Surface{sur_soil}; Delete; }{ Surface{all_sur_cable(),sur_EMdom}; };
sur_soil=news-1;

sur_airout = news; Rectangle(news) = {-dinf_th, depth_cable, 0, 2*dinf_th, dinf_th_air, 0};


all_sur() = Surface{:};
BooleanFragments{Surface{all_sur()}; Delete;}{ } // WATCH OUT: Fragments needed for forcing Coherence
bnd() = CombinedBoundary{Surface{all_sur()};};
// Printf("",bnd());

bnd_EMdom() = CombinedBoundary{Surface{sur_EMdom};};
// Printf("",bnd_EMdom());



// Adjusting some characteristic lengths
cl = dtot/5;
// Characteristic Length { Point{:} } = cl;
Characteristic Length { PointsOf{ Surface{sur_wire(),sur_semi_in(),sur_xlpe(),sur_semi_out(),sur_al()}; } } = cl/32 ;
Characteristic Length { PointsOf{ Surface{sur_ps_exact(), sur_steel_armour}; } } = cl/32 ;
Characteristic Length { PointsOf{ Surface{sur_pc, sur_steel_pipe}; } } = cl/16 ;
Characteristic Length { PointsOf{ Line{bnd_EMdom(1)};} } = 2*cl ;
Characteristic Length { PointsOf{ Surface{sur_airout}; Line{bnd(0)}; } } = 5*cl;
If(Flag_defect_in_XLPE)
  Characteristic Length { PointsOf{ Surface{sur_defect}; } } = cl/64 ;
EndIf


// ===========================================
// Physical regions => link to pro-file and FE
// ===========================================

Physical Surface("wire 1", WIRE+0) =  sur_wire(0);
Physical Surface("wire 2", WIRE+1) =  sur_wire(2);
Physical Surface("wire 3", WIRE+2) =  sur_wire(1);

Physical Surface("inner semiconductor", SEMI_IN) = sur_semi_in();
Physical Surface("XLPE", XLPE) =  sur_xlpe();
Physical Surface("outer semiconductor", SEMI_OUT) = sur_semi_out();
Physical Surface("Aluminum tape", APL) = sur_al();
Physical Surface("Air in cable",AIR_IN) = {sur_air(0), sur_air_in()};
Physical Surface("Polyethylene sheath", POLYETHYLENE_SHEATH) = sur_ps_exact();
Physical Surface("Steel armour", STEEL_ARMOUR) = sur_steel_armour;

Physical Surface("Steel pipe", STEEL_PIPE) = sur_steel_pipe;
Physical Surface("Polyethylene cover", POLYETHYLENE_COVER) = sur_pc;


Physical Surface("Soil (EM)", SOIL_EM) = sur_EMdom;
Physical Surface("Soil (thermal)", SOIL_TH) = sur_soil;
Physical Surface("Air above soil", AIR_OUT) = sur_airout;


If(Flag_defect_in_XLPE)
  Physical Surface("Defect in XLPE (phase 2)", DEFECT) = sur_defect;
  Color Cyan {Surface{sur_defect};}
EndIf


Physical Line("Outer boundary cable", OUTBND_CABLE) = 37;
Physical Line("Outer boundary (EM)", OUTBND_EM) = bnd_EMdom(1);
Physical Line("Outer boundary (TH)", OUTBND_TH) = bnd();
Physical Line("Air/Soil interface", INTERFACE_AIR_SOIL) = 159; // if force convection at air/soil interface


Color Orange {Surface{sur_soil, sur_EMdom};}

Color SteelBlue {Surface{sur_steel_armour,sur_steel_pipe};}
Color Yellow {Surface{sur_ps_exact(),sur_pc};} // polyethylene
Color Cyan {Surface{sur_air(0), sur_air_in(), sur_airout};}
Color Blue {Surface{sur_al()} ;}
Color Green {Surface{sur_xlpe()} ;}
Color Pink {Surface{sur_semi_in(),sur_semi_out()} ;}
Color Red {Surface{sur_wire()} ;}
//+
Show "*";
