//+ LowLevelFEM
//+
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {1, 0, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Transfinite Curve {1} = 2 Using Progression 1;

Mesh 1;
//+
Physical Curve("rod", 2) = {1};
//+
Physical Point("left", 3) = {1};
//+
Physical Point("right", 4) = {2};
