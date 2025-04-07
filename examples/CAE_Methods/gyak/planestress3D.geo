//+
SetFactory("OpenCASCADE");
//b=10;
//l=50;
//n=3;
Box(1) = {0, 0, 0, l, l, b};
//+
Transfinite Curve {1, 3, 7, 5} = n+1 Using Progression 1;
//+
Transfinite Curve {10, 9, 11, 12, 2, 4, 8, 6} = Ceil(l/(b/n)) Using Progression 1;
//+
Transfinite Surface {6} = {1, 5, 7, 3};
//+
Transfinite Surface {5};
//+
Transfinite Surface {1};
//+
Transfinite Surface {2};
//+
Transfinite Surface {3};
//+
Transfinite Surface {4};
//+
Transfinite Volume{1};
//+
Recombine Surface {6, 5, 4, 3, 1, 2};
//+
Mesh.ElementOrder=1;
//+
Mesh 3;
//+
Physical Volume("body", 13) = {1};
//+
Physical Surface("left", 14) = {1};
//+
Physical Surface("right", 15) = {2};
//+
Physical Surface("front", 16) = {3};
//+
Physical Surface("bottom", 17) = {5};
//+
Point(9) = {l/2, l/2, 0, 1.0};
//+
Point(10) = {l/2, l/2, b, 1.0};
//+
Line(13) = {9, 10};
//+
Physical Curve("path", 18) = {13};
//+
Point(11) = {0, l/2, b/2, 1.0};
//+
Point(12) = {l, l/2, b/2, 1.0};
//+
Line(14) = {11, 12};
//+
Physical Curve("path2", 19) = {14};
