//+
SetFactory("OpenCASCADE");
Box(1) = {0, 0, 0, 3, 30, 200};
//+
Transfinite Curve {1, 5, 3, 7} = 201 Using Progression 1;
//+
Transfinite Curve {6, 2, 4, 8} = 31 Using Progression 1;
//+
Transfinite Curve {12, 10, 9, 11} = 4 Using Progression 1;
//+
Transfinite Surface {1};
//+
Transfinite Surface {4};
//+
Transfinite Surface {2};
//+
Transfinite Surface {3};
//+
Transfinite Surface {6};
//+
Transfinite Surface {5};
//+
Transfinite Volume{1};
//+
Recombine Surface {1, 2, 4, 3, 5, 6};
//+
Mesh 3;
//+
Physical Surface("supp", 13) = {5};
//+
Physical Volume("beam", 26) = {1};
//+
Physical Surface("load", 27) = {6};

//+
Point(9) = {0, 0, 100, 1.0};
//+
Point(10) = {3, 0, 100, 1.0};
//+
Point(11) = {3, 30, 100, 1.0};
//+
Point(12) = {0, 30, 100, 1.0};
//+
Line(13) = {9, 10};
//+
Line(14) = {10, 11};
//+
Line(15) = {11, 12};
//+
Line(16) = {12, 9};
//+
Physical Curve("path", 28) = {13, 14};
//+
Physical Curve("path1", 29) = {13};
//+
Physical Curve("path2", 30) = {14};
