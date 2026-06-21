//+ LowLevelFEM
//+
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {20, 0, 0, 1.0};
//+
Point(3) = {20, 2, 0, 1.0};
//+
Point(4) = {2, 2, 0, 1.0};
//+
Point(5) = {2.2, 30, 0, 1.0};
//+
Point(6) = {0, 30, 0, 1.0};
//+
Point(7) = {0, 2, 0, 1.0};
//+
Point(8) = {-30, 2, 0, 1.0};
//+
Point(9) = {-30, -10, 0, 1.0};
//+
Point(10) = {-28, -10, 0, 1.0};
//+
Point(11) = {-28, 0, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 5};
//+
Line(5) = {5, 6};
//+
Line(6) = {6, 7};
//+
Line(7) = {7, 8};
//+
Line(8) = {8, 9};
//+
Line(9) = {9, 10};
//+
Line(10) = {10, 11};
//+
Line(11) = {11, 1};
//+
Curve Loop(1) = {6, 7, 8, 9, 10, 11, 1, 2, 3, 4, 5};
//+
Plane Surface(1) = {1};

Mesh.ElementOrder=2;
MeshSize {:} = 0.5;
Mesh 2;
//+
Physical Surface("body", 12) = {1};
//+
Physical Curve("perem", 13) = {6, 5, 4, 3, 1, 11, 10, 9, 8, 7};
//+
