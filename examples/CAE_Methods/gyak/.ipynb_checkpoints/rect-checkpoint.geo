//+
SetFactory("OpenCASCADE");
Rectangle(1) = {0, 0, 0, 100, 50, 0};
//+
Physical Curve("left", 5) = {4};
//+
Physical Curve("right", 6) = {2};
//+
Physical Surface("body", 7) = {1};

Mesh.ElementOrder = 4;
//+
MeshSize {4, 3, 2, 1} = 5;
//+
Recombine Surface {1};

Mesh 2;


