//+
SetFactory("OpenCASCADE");
//+
Box(1) = {-5, -5, 0, 10, 10, 100};
//+
MeshSize {:} = 3;
Mesh.ElementOrder=2;
Mesh 3;
//+
Physical Volume("rud", 4) = {1};
