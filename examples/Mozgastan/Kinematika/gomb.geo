//+
SetFactory("OpenCASCADE");
//+
Sphere(1) = {0, 0, 0, 1, -Pi/2, Pi/2, 2*Pi};
//+
Mesh 3;
//+
Physical Volume("volu", 4) = {1};
