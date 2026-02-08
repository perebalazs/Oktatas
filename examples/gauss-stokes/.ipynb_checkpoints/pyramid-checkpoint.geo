//+
SetFactory("OpenCASCADE");
//+
Rectangle(1) = {-0.125, -0.125, 0, 0.25, 0.25, 0};
//+
Physical Surface("square", 5) = {1};
//+
MeshSize {:} = 0.01;
Mesh.ElementOrder=4;
//+
Mesh 2;
