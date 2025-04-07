//+
SetFactory("OpenCASCADE");
//l = 100;
//d = 10;

Rectangle(1) = {-l/2, -l/2, 0, l, l, 0};
//+
Disk(2) = {0, 0, 0, d/2, d/2};
//+
BooleanDifference{ Surface{1}; Delete; }{ Surface{2}; Delete; }
//+
Physical Surface("body", 10) = {1};
//+
Physical Curve("left", 11) = {7};
//+
Physical Curve("right", 12) = {8};
//+
Physical Curve("bottom", 13) = {6};
//+
Physical Curve("circle", 14) = {5};
//+
MeshSize {5} = esh;
MeshSize {6,7,8,9} = es;
//+
Mesh 2;
//+
