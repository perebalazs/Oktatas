//+
SetFactory("OpenCASCADE");
//+
Rectangle(1) = {0, 0, 0, 1, 1, 0};
//+
Transfinite Curve {3, 2, 1, 4} = 5 Using Progression 1;
//+
Transfinite Surface {1};
//+
Recombine Surface {1};
//+
Mesh 2;
//+
Physical Surface("plate", 5) = {1};
//+
Physical Curve("supp", 6) = {4};
//+
Physical Curve("load", 7) = {2};
