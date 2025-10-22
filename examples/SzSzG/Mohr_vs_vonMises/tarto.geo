//+
SetFactory("OpenCASCADE");
Box(1) = {-30, -30, 0, 60, 60, 10};
//+
Box(2) = {-10, -10, 10, 20, 20, 110};
//+
Box(3) = {-3, -10, 120, 6, 20, 20};
//+
Coherence;
//+
Cylinder(4) = {-18, -18, -1, 0, 0, 12, 5, 2*Pi};
Cylinder(5) = { 18, -18, -1, 0, 0, 12, 5, 2*Pi};
Cylinder(6) = { 18,  18, -1, 0, 0, 12, 5, 2*Pi};
Cylinder(7) = {-18,  18, -1, 0, 0, 12, 5, 2*Pi};
//+
Cylinder(8) = {0, -11, 50, 0, 22, 0, 5, 2*Pi};
//+
Cylinder(9) = {-4, 0, 130, 8, 0, 0, 5, 2*Pi};
//+
BooleanDifference{ Volume{1}; Volume{2}; Volume{3}; Delete; }{ Volume{4}; Volume{7}; Volume{6}; Volume{5}; Volume{8}; Volume{9}; Delete; }
//+
//+
BooleanUnion{ Volume{2}; Delete; }{ Volume{3}; Delete; }
//+
BooleanUnion{ Volume{1}; Delete; }{ Volume{2}; Delete; }
//+
Fillet {1}{13, 14, 16, 15, 47, 51, 3, 19, 1, 6, 7, 2, 8, 9}{2}
//+
MeshSize {:} = 3;
Mesh.ElementOrder=2;
//+
Mesh 3;
//+
Physical Surface("base", 97) = {2};
//+
Physical Surface("hole1", 98) = {12, 14, 13, 11};
//+
Physical Surface("hole2", 99) = {39};
//+
Physical Curve("perim", 100) = {38, 34, 32, 36};
//+
Physical Volume("part", 101) = {1};
