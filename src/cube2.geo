//+
SetFactory("OpenCASCADE");
Box(1) = {0, 0, 0, 100, 10, 10};
//+
Physical Surface("supp", 13) = {1};
//+
Physical Surface("load", 14) = {2};

//MeshSize {1:4} = 5;
Mesh.ElementOrder = 4;
//Mesh.HighOrderOptimize = 2;

SetName "cube2";
Mesh 3;
Mesh.SaveAll=1;
Save "cube2.msh";

