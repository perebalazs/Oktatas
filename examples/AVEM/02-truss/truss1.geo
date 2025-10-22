//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {100, 0, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Transfinite Curve {1} = 30 Using Progression 1;
//+
Mesh 1;
//+
Physical Point("left", 2) = {1};
//+
Physical Point("right", 3) = {2};
//+
Physical Curve("rod", 4) = {1};
