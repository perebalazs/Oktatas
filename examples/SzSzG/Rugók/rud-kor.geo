//+
SetFactory("OpenCASCADE");
//+
Cylinder(1) = {0, 0, 0, 0, 0, 100, 5, 2*Pi};
//+
MeshSize {:} = 3;
Mesh.ElementOrder=2;
Mesh 3;
//+
Physical Volume("rud", 4) = {1};
