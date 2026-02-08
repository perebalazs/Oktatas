//+
SetFactory("OpenCASCADE");
//+
Torus(1) = {0, 0, 0, R, r, 2*Pi};
//+
Mesh.ElementOrder=3;
//+
MeshSize {:} = 1;
//+
Mesh 3;
//+
Physical Volume("volu", 3) = {1};
//+
Physical Surface("surf", 4) = {1};
