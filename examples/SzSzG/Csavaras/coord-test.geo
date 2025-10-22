//+
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {1, 0, 0, 1.0};
//+
Point(3) = {0, 1, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 1};
//+
Transfinite Curve {3, 2, 1} = 2 Using Progression 1;
Mesh 1;
//+
Physical Curve("truss", 4) = {1, 2, 3};
//+
Physical Point("left", 5) = {1};
//+
Physical Point("right", 6) = {2};
//+
Physical Point("top", 7) = {3};
