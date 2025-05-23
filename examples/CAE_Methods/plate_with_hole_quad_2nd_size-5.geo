//+
SetFactory("OpenCASCADE");

R = 20;

Rectangle(1) = {0, 0, 0, 100, 100, 0};
//+
Disk(2) = {0, 0, 0, R, R};
//+
BooleanDifference{ Surface{1}; Delete; }{ Surface{2}; Delete; }
//+

//+
Point(6) = {0, 0, 0, 1.0};
//+
Line(6) = {6, 4};
//+
BooleanFragments{ Surface{1}; Delete; }{ Curve{6}; Delete; }
//+
Recursive Delete {
  Curve{8}; Point{7}; 
}

MeshSize {6, 3, 5, 2, 1, 4} = 5;

Recombine Surface {2, 1};
Mesh.ElementOrder = 2;
Mesh 2;
//+
Physical Surface("body", 8) = {1, 2};
//+
Physical Curve("left", 9) = {6};
//+
Physical Curve("bottom", 10) = {4};
//+
Physical Curve("right", 11) = {3};
//+
Physical Curve("top", 12) = {7};
//+
Physical Curve("path", 13) = {2};
//+
