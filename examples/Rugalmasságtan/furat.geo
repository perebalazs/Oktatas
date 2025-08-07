//+
SetFactory("OpenCASCADE");
Rectangle(1) = {-500, -50, 0, 1000, 100, 0};
//+
Disk(2) = {0, 0, 0, 5, 5};
//+
BooleanDifference{ Surface{1}; Delete; }{ Surface{2}; Delete; }
//+
MeshSize {:} = 5;
MeshSize {5} = 0.1;
//+
Mesh.ElementOrder=4;
//+
Recombine Surface {1};
//+
Mesh 2;
//+
Physical Surface("plate", 10) = {1};
//+
Physical Curve("bottom", 11) = {6};
//+
Physical Curve("top", 12) = {9};
//+
Physical Curve("left", 13) = {7};
//+
Physical Curve("right", 14) = {8};
//+
