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
Point(1) = {0, -h/2, 0, 1.0};
//+
Point(2) = {b, -h/2, 0, 1.0};
//+
Point(3) = {b, -h/2+t, 0, 1.0};
//+
Point(4) = {t, -h/2+t, 0, 1.0};
//+
Point(5) = {t, h/2-t, 0, 1.0};
//+
Point(6) = {b, h/2-t, 0, 1.0};
//+
Point(7) = {b, h/2, 0, 1.0};
//+
Point(8) = {0, h/2, 0, 1.0};
//+
Point(9) = {0, h/2-t, 0, 1.0};
//+
Point(10) = {0, -h/2+t, 0, 1.0};
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
Line(10) = {10, 1};
//+
Line(11) = {9, 5};
//+
Line(12) = {4, 10};
//+
Curve Loop(1) = {1, 2, 3, 12, 10};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {12, -9, 11, -4};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {5, 6, 7, 8, 11};
//+
Plane Surface(3) = {3};

l1 = t / n;
n2 = Ceil((b - t) / l1);
n3 = n + n2;
n4 = Ceil((h - 2 * t) / l1);

//+
Transfinite Curve {11, 6, 12, 2, 8, 10} = n+1 Using Progression 1;
//+
Transfinite Curve {3, 5} = n2+1 Using Progression 1;
//+
Transfinite Curve {7, 1} = n3+1 Using Progression 1;
//+
Transfinite Curve {4, 9} = n4+1 Using Progression 1;
//+
Transfinite Surface {1} = {1, 2, 3, 10};
//+
Transfinite Surface {2} = {10, 4, 5, 9};
//+
Transfinite Surface {3} = {9, 6, 7, 8};
//+
Recombine Surface {1, 2, 3};

Mesh 2;

//+
Physical Surface("section", 13) = {3, 2, 1};
//+
Physical Curve("perimeter", 14) = {7, 8, 9, 10, 1, 2, 3, 4, 5, 6};
//+
Physical Point("reference", 15) = {1};
