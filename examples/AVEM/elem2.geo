//+
SetFactory("OpenCASCADE");
Rectangle(1) = {-1, -1, 0, 2, 2, 0};
//+
Transfinite Curve {3, 2, 1, 4} = 2 Using Progression 1;
//+
Transfinite Surface {1};
//+
Recombine Surface {1};

Mesh.ElementOrder=2;
Mesh 2;
