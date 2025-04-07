//+
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {100, 0, 0, 1.0};
//+
Point(3) = {100, 50, 0, 1.0};
//+
Point(4) = {55, 50, 0, 1.0};
//+
Point(5) = {50, 55, 0, 1.0};
//+
Point(6) = {55, 55, 0, 1.0};
//+
Point(7) = {50, 100, 0, 1.0};
//+
Point(8) = {0, 100, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Circle(4) = {4, 6, 5};
//+
Line(5) = {5, 7};
//+
Line(6) = {7, 8};
//+
Line(7) = {8, 1};
//+
Curve Loop(1) = {7, 1, 2, 3, 4, 5, 6};
//+
Plane Surface(1) = {1};
//+
Disk(2) = {30, 30, 0, 10, 10};
//+
BooleanDifference{ Surface{1}; Delete; }{ Surface{2}; Delete; }
//+
Transfinite Curve {8} = 16 Using Progression 1;
//+
Physical Surface("body", 8) = {1};
//+
Physical Curve("supp", 9) = {10};
//+
Physical Curve("load", 10) = {13};
//+
Mesh.ElementOrder=3;
//+
MeshSize {11, 13, 15, 12, 10} = 10;
//+
MeshSize {8, 14, 16} = 2;
//+
Mesh 2;
