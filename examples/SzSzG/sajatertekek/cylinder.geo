//+
SetFactory("OpenCASCADE");
Cylinder(1) = {0, 0, 0, 0, 0, 50, 5, 2*Pi};

MeshSize {:} = 2;
Mesh.ElementOrder=3;
//+
Mesh 3;
//+
Physical Volume("body", 4) = {1};
