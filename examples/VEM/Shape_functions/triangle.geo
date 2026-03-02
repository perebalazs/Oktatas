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
Curve Loop(1) = {3, 1, 2};
//+
Plane Surface(1) = {1};
//+
Transfinite Curve {3, 1, 2} = 2 Using Progression 1;
//+
Transfinite Surface {1};
//+
Mesh.ElementOrder=4;
//+
Mesh 2;
//+
Physical Surface("surf", 4) = {1};
