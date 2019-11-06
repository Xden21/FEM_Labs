// Gmsh project created on Wed Nov 06 14:03:06 2019
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {1, 0, 0, 1.0};
//+
Point(3) = {-1, 0, 0, 1.0};
//+
Point(4) = {0, 1, 0, 1.0};
//+
Point(5) = {0, -1, 0, 1.0};
//+
Point(6) = {1.5, 0, 0, 1.0};
//+
Point(7) = {-1.5, 0, 0, 1.0};
//+
Point(8) = {0, 1.5, 0, 1.0};
//+
Point(9) = {0, -1.5, 0, 1.0};
//+
Circle(1) = {2, 1, 4};
//+
Circle(2) = {4, 1, 3};
//+
Circle(4) = {5, 1, 2};
//+
Circle(5) = {3, 1, 5};
//+
Circle(6) = {6, 1, 8};
//+
Circle(7) = {8, 1, 7};
//+
Circle(8) = {7, 1, 9};
//+
Circle(9) = {9, 1, 6};
//+
Curve Loop(1) = {7, 8, 9, 6};
//+
Curve Loop(2) = {2, 5, 4, 1};
//+
Plane Surface(1) = {1, 2};
//+
Extrude {0, 0, 10} {
  Surface{1};
}
