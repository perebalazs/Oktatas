//+
SetFactory("OpenCASCADE");

//+
Point(1) = {-1, -0, 0, 1.0};
//+
Point(2) = {-1, -1, 0, 1.0};
//+
Point(3) = {-0, -0, 0, 1.0};
//+
Point(4) = {1, 1, 0, 1.0};
//+
Point(5) = {-1, 1, 0, 1.0};
//+
Point(6) = {1, -1, 0, 1.0};
//+
Point(7) = {1, -0, 0, 1.0};
//+
Point(8) = {-0, 1, 0, 1.0};
//+
Point(9) = {0, -1, 0, 1.0};
//+
Line(1) = {2, 9};
//+
Line(2) = {9, 6};
//+
Line(3) = {6, 7};
//+
Line(4) = {7, 4};
//+
Line(5) = {4, 8};
//+
Line(6) = {8, 5};
//+
Line(7) = {5, 1};
//+
Line(8) = {1, 2};
//+
Line(9) = {9, 3};
//+
Line(10) = {3, 8};
//+
Line(11) = {1, 3};
//+
Line(12) = {3, 7};


//+
Mesh 3;

//+
Physical Point("A", 2) = {1};
//+
Physical Point("B", 3) = {7};
//+
Physical Point("C", 4) = {8};
//+
Physical Curve("frame", 13) = {2, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 1};
//+
Physical Point("O", 14) = {3};
//+
Physical Point("D", 15) = {4};
