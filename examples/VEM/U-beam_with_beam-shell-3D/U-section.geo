//+ LowLevelFEM
//+
SetFactory("OpenCASCADE");
//+
b = DefineNumber[ 50, Name "Parameters/base" ];
//+
h = DefineNumber[ 100, Name "Parameters/height" ];
//+
t = DefineNumber[ 5, Name "Parameters/thikness" ];
//+
n = DefineNumber[ 3, Name "Parameters/elements_per_edge" ];
//+
Point(1) = {-b/2, -h/2, 0, 1.0};
//+
Point(2) = {b/2, -h/2, 0, 1.0};
//+
Point(3) = {b/2, -h/2+t, 0, 1.0};
//+
Point(4) = {t/2, -h/2+t, 0, 1.0};
//+
Point(5) = {t/2, h/2-t, 0, 1.0};
//+
Point(6) = {b/2, h/2-t, 0, 1.0};
//+
Point(7) = {b/2, h/2, 0, 1.0};
//+
Point(8) = {-b/2, h/2, 0, 1.0};
//+
Point(9) = {-b/2, h/2-t, 0, 1.0};
//+
Point(10) = {-t/2, h/2-t, 0, 1.0};
//+
Point(11) = {-t/2, -h/2+t, 0, 1.0};
//+
Point(12) = {-b/2, -h/2+t, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 5};
//+
Line(5) = {5, 6};
//+
Line(6) = {6, 7};
//+
Line(7) = {7, 8};
//+
Line(8) = {8, 9};
//+
Line(9) = {9, 10};
//+
Line(10) = {10, 11};
//+
Line(11) = {11, 12};
//+
Line(12) = {12, 1};
//+
Curve Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
//+
Plane Surface(1) = {1};

n1 = n;
l1 = t / n;
n2 = Ceil((b - t) / 2 / l1) + 1;
n3 = 2 * n2 + n;
n4 = Ceil((h - 2 * t) / l1) + 1;
//+
Transfinite Curve {8, 6, 2, 12} = n + 1 Using Progression 1;
//+
Transfinite Curve {9, 5, 11, 3} = n2 + 1 Using Progression 1;
//+
Transfinite Curve {7, 1} = n3 + 1 Using Progression 1;
//+
Transfinite Curve {4, 10} = n4 Using Progression 1;


//+
//Transfinite Surface {1};
//+
//Recombine Surface {1};

Mesh 2;



//+
Physical Surface("section", 13) = {1};
//+
Physical Curve("perimeter", 14) = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
//+
Physical Point("reference", 15) = {1};
