//+
SetFactory("OpenCASCADE");
//+
Circle(1) = {0, 0, 0, 25, 0, Pi/2};
//+
Circle(2) = {0, 0, 0, 23, 0, Pi/2};
//+
Line(3) = {4, 2};
//+
Line(4) = {3, 1};
//+
Curve Loop(1) = {2, 3, -1, -4};
//+
Plane Surface(1) = {1};
//+
Extrude {0, 0, 100} {
  Surface{1}; 
}
//+
Extrude {0, 0, 100} {
  Surface{6}; 
}
//+
Rotate {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Duplicata { Volume{1,2}; }
}
//+
Rotate {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Duplicata { Volume{3,4}; }
}
//+
Rotate {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Duplicata { Volume{5,6}; }
}
Coherence;
//+
Box(9) = {-0.001, -20, -1, 0.001, -10, 201};
//+
BooleanDifference{ Volume{3}; Volume{4}; Volume{2}; Volume{1}; Volume{6}; Volume{5}; Volume{7}; Volume{8}; Delete; }{ Volume{9}; Delete; }
//+
Transfinite Curve {17, 7, 6, 32, 26, 25, 20, 12, 11, 45, 58, 39, 57, 38, 64} = 4 Using Progression 1;
//+
Transfinite Curve {16, 14, 13, 18, 29, 31, 42, 54, 44, 59, 5, 2, 1, 8, 21, 24, 34, 62, 37, 65} = 20 Using Progression 1;
//+
Transfinite Curve {15, 19, 4, 10, 3, 9, 30, 33, 23, 28, 22, 27, 55, 60, 53, 56, 61, 63, 51, 52, 48, 50, 47, 49} = 24 Using Progression 1;
//+
Transfinite Surface {:};
//+
Transfinite Volume{:};
//+
Recombine Surface {:};
ReverseMesh Surface{38};
Mesh 3;
//+//+
Physical Volume("body", 66) = {4, 6, 2, 8, 3, 5, 1, 7};
//+
Physical Surface("left", 67) = {15, 43, 32, 5};
//+
Physical Surface("right", 68) = {11, 20, 40, 36};
//+
Physical Surface("mid", 69) = {16, 38, 33, 6};
//+
Physical Curve("path", 70) = {12};
