//+
SetFactory("OpenCASCADE");
Box(1) = {0, 0, 0, 10, 10, 10};
//+
Recombine Surface {2, 1, 4, 5, 6, 3};
//+
Transfinite Curve {9, 8, 5, 11, 4, 1, 10, 6, 7, 3, 12, 2} = 15 Using Progression 1;
//+
Transfinite Surface {1:6};
//+
Transfinite Surface {3};
//+
Transfinite Surface {5};
//+
Transfinite Volume{1};

Mesh 3;

//+
Physical Volume("cube", 13) = {1};
//+
Physical Surface("supp", 14) = {5};
//+
Physical Surface("load", 15) = {6};
//+
Point(9) = {5, 0, 5, 1.0};
//+
Point(10) = {5, 10, 5, 1.0};
//+
Line(13) = {9, 10};
//+
Physical Curve("path", 16) = {13};
