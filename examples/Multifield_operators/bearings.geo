//+
SetFactory("OpenCASCADE");

//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {100, 0, 0, 1.0};
//+
Point(3) = {100, 30, 0, 1.0};
//+
Point(4) = {0, 30, 0, 1.0};
//+
Point(5) = {5, 4, 0, 1.0};
//+
Point(6) = {95, 4, 0, 1.0};
//+
Point(7) = {95, 6, 0, 1.0};
//+
Point(8) = {5, 6, 0, 1.0};
//+
Point(9) = {5, 14, 0, 1.0};
//+
Point(10) = {95, 14, 0, 1.0};
//+
Point(11) = {95, 16, 0, 1.0};
//+
Point(12) = {5, 16, 0, 1.0};
//+
Point(13) = {5, 24, 0, 1.0};
//+
Point(14) = {95, 24, 0, 1.0};
//+
Point(15) = {95, 26, 0, 1.0};
//+
Point(16) = {5, 26, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Line(5) = {5, 6};
//+
Line(6) = {6, 7};
//+
Line(7) = {7, 8};
//+
Line(8) = {8, 5};
//+
Line(9) = {9, 10};
//+
Line(10) = {10, 11};
//+
Line(11) = {11, 12};
//+
Line(12) = {12, 9};
//+
Line(13) = {13, 14};
//+
Line(14) = {14, 15};
//+
Line(15) = {15, 16};
//+
Line(16) = {16, 13};
//+
Curve Loop(1) = {4, 1, 2, 3};
//+
Curve Loop(2) = {5, 6, 7, 8};
//+
Curve Loop(3) = {9, 10, 11, 12};
//+
Curve Loop(4) = {13, 14, 15, 16};
//+
Plane Surface(1) = {1, 2, 3, 4};
//+
Curve Loop(5) = {5, 6, 7, 8};
//+
Plane Surface(2) = {5};
//+
Curve Loop(6) = {9, 10, 11, 12};
//+
Plane Surface(3) = {6};
//+
Curve Loop(7) = {13, 14, 15, 16};
//+
Plane Surface(4) = {7};

Mesh 2;

//+
Physical Surface("rubber", 17) = {1};
//+
Physical Surface("steel", 18) = {4, 3, 2};
//+
Physical Curve("left", 19) = {4};
//+
Physical Curve("right", 20) = {2};
//+
Physical Curve("bottom", 21) = {1};
//+
Physical Curve("top", 22) = {3};
