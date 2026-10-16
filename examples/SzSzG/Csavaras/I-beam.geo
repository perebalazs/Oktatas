//+ LowLevelFEM
//+
SetFactory("OpenCASCADE");
//+
b = DefineNumber[ 100, Name "Parameters/base" ];
//+
h = DefineNumber[ 200, Name "Parameters/hight" ];
//+
t = DefineNumber[ 10, Name "Parameters/tickness" ];
//+
L = DefineNumber[ 300, Name "Parameters/length" ];
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
Line(13) = {11, 4};
//+
Line(14) = {10, 5};
//+
Curve Loop(1) = {7, 8, 9, 14, 5, 6};
//+
Surface(1) = {1};
//+
Curve Loop(3) = {4, -14, 10, 13};
//+
Surface(2) = {3};
//+
Curve Loop(5) = {1, 2, 3, -13, 11, 12};
//+
Surface(3) = {5};
//+
n = DefineNumber[ 3, Name "Parameters/n" ];
es = t / n;
//+
Transfinite Curve {8, 6, 14, 13, 12, 2} = Ceil(t/es)+1 Using Progression 1;
//+
Transfinite Curve {10, 4} = Ceil((h-2*t)/es)+1 Using Progression 1;
//+
Transfinite Curve {9, 5, 11, 3} = Ceil(((b-t)/2/es)+1) Using Progression 1;
//+
Transfinite Curve {7, 1} = 2*Ceil(((b-t)/2/es)+1)+n-1 Using Progression 1;
//+
Transfinite Surface {3} = {1, 2, 3, 12};
//+
Transfinite Surface {1} = {9, 6, 7, 8};
//+
Transfinite Surface {2} = {11, 4, 5, 10};
//+
Recombine Surface {3, 2, 1};
//+
Extrude {0, 0, L/2} {
  Surface{1}; Surface{2}; Surface{3}; Layers {Ceil(L/2 / (t/n))}; Recombine;
}
//+
Extrude {0, 0, L/2} {
  Surface{10}; Surface{14}; Surface{20}; Layers {Ceil(L/2 / (t/n))}; Recombine;
}
Periodic Surface{1}={27} Translate {0, 0, -L};
Periodic Surface{2}={31} Translate {0, 0, -L};
Periodic Surface{3}={37} Translate {0, 0, -L};
//+
Point(37) = {0, 0, L, 1.0};
//+
p = DefineNumber[ 1, Name "Parameters/order" ];
Mesh.ElementOrder=p;
//+
Mesh 3;

//+
Physical Volume("beam", 67) = {1, 2, 3, 4, 5, 6};
//+
Physical Volume("leftbeam", 68) = {1, 2, 3};
//+
Physical Surface("left", 69) = {1, 2, 3};
//+
Physical Surface("right", 70) = {27, 31, 37};
//+
Physical Point("remote", 71) = {37};
//+
Physical Surface("mid", 72) = {10, 14, 20};
//+
Physical Point("A", 73) = {1};
