//+
Point(1) = {0, 0, 0, 20};
//+
Point(2) = {10, 0, 0, 20};
//+
Point(3) = {20, 0, 0, 20};
//+
Point(4) = {10, 10, 0, 20};
//+
Point(5) = {0, 10, 0, 20};
//+
Line(1) = {5, 1};
//+
Line(2) = {1, 2};
//+
Line(3) = {2, 3};
//+
Line(4) = {3, 4};
//+
Line(5) = {4, 5};
//+
Line(6) = {2, 4};

//+
Curve Loop(1) = {1, 2, 6, 5};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {4, -6, 3};
//+
Plane Surface(2) = {2};
//
Mesh.ElementOrder = 1;
Recombine Surface{1};
Mesh 2;

View "field" {SQ (0,0,0, 10,0,0, 10,10,0, 0,10,0) {1,1,1,1};};
View "field" {ST (10,0,0, 20,0,0, 10,10,0) {2,2,2};};
//+
//+
Physical Curve("supp", 7) = {1};
//+
Physical Curve("load", 8) = {4};
Mesh.SaveAll=1;
Save "postproc-test.msh";

