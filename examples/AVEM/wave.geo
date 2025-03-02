//+
SetFactory("OpenCASCADE");
Rectangle(1) = {0, 0, 0, 100, 10, 0};
//+
Transfinite Curve {4, 2} = 11 Using Progression 1;
//+
Transfinite Curve {3, 1} = 101 Using Progression 1;
//+
Transfinite Surface {1};
//+
Recombine Surface {1};
//+
Mesh.ElementOrder = 2;
//+
Mesh 2;
//+
Physical Curve("supp", 5) = {2};
//+
Physical Surface("body", 6) = {1};
