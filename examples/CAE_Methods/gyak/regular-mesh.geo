//+ LowLevelFEM
//+
SetFactory("OpenCASCADE");
//+
Rectangle(1) = {0, 0, 0, 200, 100, 0};
//+
MeshSize {1, 4} = 5;
//+
MeshSize {3, 2} = 20;
//+
Transfinite Curve {4, -2} = 11 Using Progression 1.1;
//+
Transfinite Curve {3, 1} = 21 Using Progression 1;
//+
Transfinite Surface {1} = {1, 2, 3, 4};
//+
Recombine Surface {1};

Mesh 2;
//+
Point(5) = {100, -10, 0, 1.0};
//+
Point(6) = {100, 110, 0, 1.0};
//+
Line(5) = {5, 6};
//+
BooleanFragments{ Surface{1}; Delete; }{ Curve{5}; Delete; }

Transfinite Curve {-2, 4,7,5} = 11 Using Progression 1.1;
//+
Transfinite Curve {1,3,6} = 21 Using Progression 1;
//+
Transfinite Surface {1} = {2,4,3,1};
//+
Recombine Surface {1};

Mesh 2;