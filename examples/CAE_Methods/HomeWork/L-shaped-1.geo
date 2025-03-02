//R=5;
//+
//es=10;
//+
Point(1) = {0, 0, 0, es};
//+
Point(2) = {100, 0, 0, es};
//+
Point(3) = {100, 50, 0, es};
//+
Point(4) = {50+R, 50, 0, esR};
//+
Point(5) = {50, 50+R, 0, esR};
//+
Point(6) = {50, 100, 0, es};
//+
Point(7) = {0, 100, 0, es};
//+
Point(8) = {50+R, 50+R, 0, 1};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Circle(4) = {4, 8, 5};
//+
Line(5) = {5, 6};
//+
Line(6) = {6, 7};
//+
Line(7) = {7, 1};
//+
Curve Loop(1) = {1, 2, 3, 4, 5, 6, 7};
//+
Plane Surface(1) = {1};
//+
Mesh 2;
//+
Physical Surface("body", 8) = {1};
//+
Physical Curve("support", 9) = {6};
//+
Physical Curve("load", 10) = {2};
//+
Physical Curve("path1", 11) = {3, 4, 5};
//+
Line(8) = {1, 8};
//+
Physical Curve("path2", 12) = {8};
