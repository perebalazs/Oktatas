//+
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {100, 0, 0, 1.0};
//+
Point(3) = {100, 8, 0, 1.0};
//+
Point(4) = {8, 8, 0, 1.0};
//+
Point(5) = {8, 192, 0, 1.0};
//+
Point(6) = {100, 192, 0, 1.0};
//+
Point(7) = {100, 200, 0, 1.0};
//+
Point(8) = {0, 200, 0, 1.0};
//+
Point(9) = {8, 0, 0, 1.0};
//+
Point(10) = {8, 200, 0, 1.0};
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
Line(8) = {8, 1};
//+
Line(9) = {9, 4};
//+
Line(10) = {10, 5};
//+
Curve Loop(1) = {7, 8, 1, 2, 3, 4, 5, 6};
//+
Plane Surface(1) = {1};
//+
BooleanFragments{ Surface{1}; Delete; }{ Curve{10}; Curve{9}; Delete; }
//+
Transfinite Curve {16} = 101 Using Progression 1;
//+
Transfinite Curve {14} = 93 Using Progression 1;
//+
Transfinite Curve {10, 9, 15, 17, 11, 19} = 5 Using Progression 1;
//+
Transfinite Curve {18, 20, 12, 13} = 47 Using Progression 1;
//+
Transfinite Surface {3};
//+
Transfinite Surface {1};
//+
Transfinite Surface {2} = {6, 7, 4, 8};
//+
Recombine Surface {1, 2, 3};
//+
Coherence;
//+
Mesh 2;
//+
Extrude {0, 0, 200} {
  Surface{1}; Surface{2}; Surface{3}; Layers {100}; Recombine;
}
//+
Mesh 3;
//+
Physical Surface("supp", 43) = {1, 2, 3};
//+
Physical Surface("load", 44) = {8, 14, 18};
//+
Physical Volume("body", 45) = {1, 2, 3};
//+
Point(21) = {100, 0, 100, 1.0};
//+
Point(22) = {0, 0, 100, 1.0};
//+
Point(23) = {0, 200, 100, 1.0};
//+
Point(24) = {100, 200, 100, 1.0};
//+
Point(25) = {100, 4, 100, 1.0};
//+
Point(26) = {4, 4, 100, 1.0};
//+
Point(27) = {4, 196, 100, 1.0};
//+
Point(28) = {100, 196, 100, 1.0};
//+
Point(29) = {100, 8, 100, 1.0};
//+
Point(30) = {8, 8, 100, 1.0};
//+
Point(31) = {8, 192, 100, 1.0};
//+
Point(32) = {100, 192, 100, 1.0};
//+
Point(33) = {0, 50, 100, 1.0};
//+
Point(34) = {0, 100, 100, 1.0};
//+
Point(35) = {0, 150, 100, 1.0};
//+
Point(36) = {8, 50, 100, 1.0};
//+
Point(37) = {8, 100, 100, 1.0};
//+
Point(38) = {8, 150, 100, 1.0};
//+
Point(39) = {50, 0, 100, 1.0};
//+
Point(40) = {50, 8, 100, 1.0};
//+
Point(41) = {50, 192, 100, 1.0};
//+
Point(42) = {50, 200, 100, 1.0};
//+
Line(43) = {21, 22};
//+
Line(44) = {22, 23};
//+
Line(45) = {23, 24};
//+
Line(46) = {25, 26};
//+
Line(47) = {26, 27};
//+
Line(48) = {27, 28};
//+
Line(49) = {29, 30};
//+
Line(50) = {30, 31};
//+
Line(51) = {31, 32};
//+
Line(52) = {33, 36};
//+
Line(53) = {34, 37};
//+
Line(54) = {35, 38};
//+
Line(55) = {39, 40};
//+
Line(56) = {41, 42};
//+
Physical Curve("path1", 57) = {43, 44, 45};
//+
Physical Curve("path2", 58) = {49, 50, 51};
//+
Physical Curve("path3", 59) = {46, 47, 48};
//+
Physical Curve("path4", 60) = {52};
//+
Physical Curve("path5", 61) = {53};
//+
Physical Curve("path6", 62) = {54};
//+
Physical Curve("path7", 63) = {55};
//+
Physical Curve("path8", 64) = {56};
