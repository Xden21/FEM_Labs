mm = 1e-3; // Unit

// Mesh sizes todo see example data.geo => they use other definitions
/* mxs = 0.05;
ms = 0.1;
mm = 1;
ml = 10; */

ppm = '{1Geometrical data/11mesh control/0';
colorro  = "LightGrey";

// Makes it possible to see inside onelab
DefineConstant[
  md = { 1.,  Name StrCat[ppm, "0Mesh density"]},
  nn_wcore   = { Ceil[md*4],  Name StrCat[ppm, "0core width"], ReadOnly 1, Highlight Str[colorro], CLose 1},
  nn_airgap  = { Ceil[md*1], Name StrCat[ppm, "1air gap width"], ReadOnly 1, Highlight Str[colorro]},
  nn_ri = { Ceil[md*12], Name StrCat[ppm,"2"], Label "1/4 shell in", ReadOnly 1, Highlight Str[colorro]},
  nn_ro = { Ceil[md*12], Name StrCat[ppm,"3"], Label "1/4 shell out", ReadOnly 1, Highlight Str[colorro]}
];

// local constants
// IE-Core measurements
wI = 24.5*mm; // Width core

wE = wI/8; // Width lateral legs
//wEc = wI/4; // Width Central legs => not necessary since we only draw half of inductor
wC = wI/4; // Width coils

hI = wE; // Height return path
hE = 5*wI/8; // Height E core
hC = wI/2; // Height coils
ag = 0.3175*mm; // length air gap

htot = hI + hE + ag; // Total height of inductor

// shell radius
Rint = 30*mm; // Internal radius
Rext = 40*mm; // External radius
