es = 20;
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {20, 0, 0, 1.0};
//+
Point(3) = {20, 10, 0, 1.0};
//+
Point(4) = {18, 10, 0, 1.0};
//+
Point(5) = {0, 10, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 5};
//+
Line(5) = {5, 1};
//+
Curve Loop(1) = {4, 5, 1, 2, 3};
//+
Plane Surface(1) = {1};
//+
Transfinite Surface {1} = {1, 2, 3, 5};
//+
Recombine Surface {1};
//+
Mesh.ElementOrder=3;
//+
Mesh 2;
//+
Physical Surface("plate", 6) = {1};
//+
Physical Curve("supp", 7) = {5};
//+
Physical Curve("load", 8) = {3};
