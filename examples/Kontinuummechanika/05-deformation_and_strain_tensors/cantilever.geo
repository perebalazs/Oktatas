//+
SetFactory("OpenCASCADE");
//+
Box(1) = {-5, -5, 0, 10, 10, 100};
//+
Transfinite Curve {12, 2, 10, 6, 11, 4, 9, 8} = 6 Using Progression 1;
//+
Transfinite Curve {3, 1, 5, 7} = 51 Using Progression 1;
//+
Transfinite Surface {:};
//+
Transfinite Volume{1};
//+
Recombine Surface {:};
//+
Mesh 3;
//+
Physical Volume("body", 13) = {1};
