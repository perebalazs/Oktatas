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


Transfinite Curve {30, 33, 23, 28, 22, 27, 15, 19, 4, 10, 3, 9, 51, 52, 48, 50, 47, 49, 43, 46, 36, 41, 35, 40} = 24 Using Progression 1;
//+
Transfinite Curve {32, 26, 25, 20, 12, 11, 17, 7, 6, 45, 39, 38} = 4 Using Progression 1;
//+
Transfinite Curve {31, 29, 14, 16, 42, 44, 13, 18, 24, 21, 2, 5, 34, 37, 1, 8} = 20 Using Progression 1;

//+
Transfinite Surface {:};
//+
Transfinite Volume{:};
//+
Recombine Surface {:};

Mesh 3;
//+//+
Physical Volume("body", 53) = {3, 4, 6, 2, 8, 5, 1, 7};
//+
Physical Surface("left", 54) = {15, 5, 24, 32};
//+
Physical Surface("right", 55) = {11, 36, 29, 20};
//+
Physical Surface("mid", 56) = {16, 6, 33, 25};
//+
Physical Curve("path", 57) = {12};
