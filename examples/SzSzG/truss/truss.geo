//+ LowLevelFEM
//+
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {100, 0, 0, 1.0};
//+
Point(3) = {50, Sqrt(3) * 50, 0, 1.0};
//+
Point(4) = {150, Sqrt(3) * 50, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 4};
//+
Line(3) = {4, 3};
//+
Line(4) = {3, 1};
//+
Line(5) = {3, 2};
//+
Transfinite Curve {:} = 2 Using Progression 1;
//+
Mesh 1;
//+
Physical Curve("truss", 6) = {1, 2, 3, 4, 5};
//+
Physical Point("A", 7) = {1};
//+
Physical Point("B", 8) = {2};
//+
Physical Point("C", 9) = {3};
//+
Physical Point("D", 10) = {4};
