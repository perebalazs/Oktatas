//+ LowLevelFEM
//+
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 10.0};
//+
Point(2) = {100, 0, 0, 10.0};
//+
Point(3) = {200, 0, 0, 10.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Mesh.ElementOrder=2;
Mesh 1;

Physical Curve("body", 3) = {1, 2};
//+
Physical Point("left", 4) = {1};
//+
Physical Point("right", 5) = {3};
//+
Physical Point("mid", 6) = {2};
