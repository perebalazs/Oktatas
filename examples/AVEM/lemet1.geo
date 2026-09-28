//+ LowLevelFEM
//+
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {100, 0, 0, 1.0};
//+
Point(3) = {100, 50, 0, 1.0};
//+
Point(4) = {0, 50, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Circle(5) = {50, 25, 0, 10, 0, 2*Pi};
//+
Curve Loop(1) = {3, 4, 1, 2};
//+
Curve Loop(2) = {5};
//+
Plane Surface(1) = {1, 2};

Mesh 2;
//+
Physical Curve("left", 6) = {4};
//+
Physical Curve("right", 7) = {2};
//+
Physical Surface("plate", 8) = {1};
