//+
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {400, 0, 0, 1.0};
//+
Point(3) = {800, 0, 0, 1.0};
//+
Point(4) = {800, 300, 0, 1.0};
//+
Point(5) = {0, 300, 0, 1.0};
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
Curve Loop(1) = {5, 1, 2, 3, 4};
//+
Plane Surface(1) = {1};
//+
Transfinite Surface {1} = {1, 3, 4, 5};
//+
MeshSize {5, 1, 2, 3, 4} = 10;
//+
Recombine Surface {1};
//+
Mesh 2;
//+

//+
Physical Surface("plate", 6) = {1};
//+
Physical Point("left", 7) = {1};
//+
Physical Point("right", 8) = {3};
//+
Physical Point("load", 9) = {2};
