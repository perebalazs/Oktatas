//+
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {0, -1000, 0, 1.0};
//+
Point(3) = {2000, 0, 0, 1.0};
//+
Line(1) = {1, 3};
//+
Line(2) = {2, 3};
//+
Transfinite Curve {1, 2} = 2 Using Progression 1;
//+
Mesh 1;
//+
Physical Curve("L1", 3) = {1};
//+
Physical Curve("L2", 4) = {2};
//+
Physical Point("P1", 5) = {1};
//+
Physical Point("P2", 6) = {2};
//+
Physical Point("P3", 7) = {3};
//+
Physical Curve("all", 8) = {1,2};
