//+
SetFactory("OpenCASCADE");
//+
Point(1) = {-50, -100, 0, 1.0};
//+
Point(2) = {50, -100, 0, 1.0};
//+
Point(3) = {50, -90, 0, 1.0};
//+
Point(4) = {5, -90, 0, 1.0};
//+
Point(5) = {5, 90, 0, 1.0};
//+
Point(6) = {50, 90, 0, 1.0};
//+
Point(7) = {50, 100, 0, 1.0};
//+
Point(8) = {-50, 100, 0, 1.0};
//+
Point(9) = {-50, 90, 0, 1.0};
//+
Point(10) = {-5, 90, 0, 1.0};
//+
Point(11) = {-5, -90, 0, 1.0};
//+
Point(12) = {-50, -90, 0, 1.0};
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
//+
Line(13) = {10, 5};
//+
Line(14) = {11, 4};
//+
Curve Loop(1) = {7, 8, 9, 13, 5, 6};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {10, 14, 4, -13};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {11, 12, 1, 2, 3, -14};
//+
Plane Surface(3) = {3};
//+
MeshSize {8, 9, 10, 5, 7, 6, 12, 1, 11, 4, 3, 2} = 3.3;
//+
Recombine Surface {1, 2, 3};
//+
Transfinite Surface {1} = {9, 6, 7, 8};
//+
Transfinite Surface {2} = {11, 4, 5, 10};
//+
Transfinite Surface {3} = {1, 2, 3, 12};
//+
//Mesh 2;
//+
Extrude {0, 0, 500} {
  Surface{3}; Surface{2}; Surface{1}; Layers {166}; Recombine;
}
//+
Extrude {0, 0, 500} {
  Surface{20}; Surface{14}; Surface{10}; Layers {166}; Recombine;
}
Coherence;
Mesh 3;
//+
Physical Volume("body", 67) = {3, 2, 1, 6, 5, 4};
//+
Physical Surface("left", 68) = {3, 2, 1};
//+
Physical Surface("right", 69) = {37, 31, 27};
//+
Physical Surface("mid", 70) = {10, 14, 20};
//+
Point(37) = {50, 95, 500, 1.0};
//+
Point(38) = {0, 95, 500, 1.0};
//+
Point(39) = {0, -95, 500, 1.0};
//+
Point(40) = {50, -95, 500, 1.0};
//+
Line(67) = {37, 38};
//+
Line(68) = {38, 39};
//+
Line(69) = {39, 40};
//+
Physical Curve("path", 71) = {67, 68, 69};
