//+
SetFactory("OpenCASCADE");
//+
Rectangle(1) = {0, 0, 0, 100, 10, 0};
//+
Physical Curve("supp", 5) = {4};
//+
Physical Curve("load", 6) = {2};
//+
Physical Surface("body", 7) = {1};

MeshSize {1:4} = 3;
Mesh.ElementOrder = 1;
//Mesh.HighOrderOptimize = 2;

SetName "bending2D";
Mesh 2;
// Mesh.SaveAll=1;
// Save "bending3D.msh";
