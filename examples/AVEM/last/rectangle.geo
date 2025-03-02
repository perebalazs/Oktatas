//+
SetFactory("OpenCASCADE");
Rectangle(1) = {0, -5, 0, 100, 10, 0};
//+
Point(5) = {0, 0, 0, 1.0};
//+
Point(6) = {100, 0, 0, 1.0};
//+
Line(5) = {5, 6};
//+
BooleanFragments{ Surface{1}; Delete; }{ Curve{5}; Delete; }
//+
Transfinite Curve {10, 5, 8} = 101 Using Progression 1;
//+
Transfinite Curve {9, 6, 11, 7} = 6 Using Progression 1;
//+
Transfinite Surface {2};
//+
Transfinite Surface {1};
//+
Recombine Surface {2, 1};
//+
Mesh 2;
//+
Physical Surface("body", 12) = {2, 1};
//+
Physical Curve("left", 13) = {6, 9};
//+
Physical Curve("right", 14) = {11, 7};
//+
Physical Point("L", 15) = {1};
//+
Point(7) = {70, -5, 0, 1.0};
//+
Point(8) = {70, 5, 0, 1.0};
//+
Line(12) = {7, 8};
//+
Physical Curve("path", 16) = {12};
