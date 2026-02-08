//+
SetFactory("OpenCASCADE");
//+
R=10;
//+
Point(1) = {R, 0, 0, 1.0};
//+
Point(2) = {0, R, 0, 1.0};
//+
Point(3) = {0, 0, R, 1.0};
//+
Point(4) = {0, 0, 0, 1.0};
//+
Circle(1) = {1, 4, 2};
//+
Circle(2) = {2, 4, 3};
//+
Circle(3) = {3, 4, 1};
//+
Curve Loop(1) = {2, 3, 1};
//+
Surface(1) = {1};
//+
Line(4) = {4, 1};
//+
Line(5) = {4, 2};
//+
Line(6) = {4, 3};
//+
Curve Loop(3) = {5, 2, -6};
//+
Plane Surface(2) = {3};
//+
Curve Loop(4) = {6, 3, -4};
//+
Plane Surface(3) = {4};
//+
Curve Loop(5) = {4, 1, -5};
//+
Plane Surface(4) = {5};
//+
Surface Loop(1) = {2, 4, 3, 1};
//+
Volume(1) = {1};


MeshSize {:} = 1;
Mesh.ElementOrder=2;
Mesh 3;

//+
Physical Curve("peri", 4) = {1, 2, 3};
//+
Physical Surface("surf", 5) = {1};
//+
Physical Volume("volu", 7) = {1};
//+
Physical Curve("line", 8) = {6, 4, 5};
