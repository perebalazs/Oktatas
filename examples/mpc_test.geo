//+
SetFactory("OpenCASCADE");
Rectangle(1) = {0, 0, 0, 50, 10, 0};
//+
Transfinite Curve {4, 2} = 5 Using Progression 1;
//+
Transfinite Curve {3, 1} = 21 Using Progression 1;
//+
Transfinite Surface {1};
//+
Recombine Surface {1};

//+
Physical Surface("body", 5) = {1};
//+
Physical Curve("left", 6) = {4};
//+
Physical Curve("right", 7) = {2};

Mesh 2;
