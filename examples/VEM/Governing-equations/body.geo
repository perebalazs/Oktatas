//+
SetFactory("OpenCASCADE");
//+
Box(1) = {0, 0, 0, 10, 10, 30};
//+
MeshSize {:} = 2;
//+
Mesh 3;
//+
Physical Volume("body", 13) = {1};
//+
Physical Surface("supp", 14) = {5};
//+
Physical Surface("load1", 15) = {6};
//+
Physical Surface("load2", 16) = {4};
