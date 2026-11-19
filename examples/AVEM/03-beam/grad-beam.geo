//+ LowLevelFEM
//+
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {1, 1, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Mesh 1;

Physical Curve("line", 2) = {1};
