//+
SetFactory("OpenCASCADE");
Box(1) = {0, 0, 0, 30, 3, 3};
//+
Transfinite Curve {1, 5, 3, 7} = 4 Using Progression 1;
//+
Transfinite Curve {6, 2, 4, 8} = 4 Using Progression 1;
//+
Transfinite Curve {12, 10, 9, 11} = 31 Using Progression 1;
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
Physical Surface("supp", 13) = {1};
//+
Physical Surface("load", 14) = {2};
//+
Physical Volume("beam", 26) = {1};

//+
Point(9) = {0, 1.5, 1.5, 1.0};
//+
Point(10) = {30, 1.5, 1.5, 1.0};
//+
Line(13) = {9, 10};
//+
Physical Curve("centerline", 27) = {13};
