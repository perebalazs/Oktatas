//+
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
Point(2) = {100, 0, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Transfinite Curve {1} = 3 Using Progression 1;

Mesh 1;
//+
Physical Point("supp", 2) = {1};
//+
Physical Point("load", 3) = {2};
//+
Physical Curve("rod", 4) = {1};
