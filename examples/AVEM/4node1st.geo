//+
Point(1) = {-0.5, 0, 0, 1.0};
//+
Point(2) = {0.5, 0, 0, 1.0};
//+
Point(3) = {0, 0.866, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 1};
//+
Curve Loop(1) = {3, 1, 2};
//+
Plane Surface(1) = {1};

Mesh.ElementOrder= 2;
Mesh 2;
//+
Physical Surface("body", 4) = {1};
