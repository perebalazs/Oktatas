//+
SetFactory("OpenCASCADE");
Box(1) = {0, 0, 0, 100, 10, 10};
//Box(2) = {20, 0, 0, 10, 10, 10};
//+
Physical Surface("supp", 13) = {1};
//+
Physical Surface("load", 14) = {2};
//+
Physical Volume("body", 15) = {1};

//+
//MeshSize {3, 7, 8, 4, 5, 6, 2, 1} = 10;
MeshSize {1:8} = 10;
Mesh.ElementOrder = 3;
//Mesh.HighOrderOptimize = 2;

SetName "bending3D";
Mesh 3;
// Mesh.SaveAll=1;
// Save "bending3D.msh";

//+
Point(9) = {10, 0, 5, 1.0};
//+
Point(10) = {10, 10, 5, 1.0};
//+
Line(13) = {9, 10};
//+
Physical Curve("path", 16) = {13};
