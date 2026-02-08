//+
SetFactory("OpenCASCADE");
//+
Box(1) = {0, 0, 0, 2, 2, 2};
//+
Mesh 3;
//+
Physical Volume("volu", 13) = {1};
//+
Physical Surface("+xy", 14) = {6};
//+
Physical Surface("-xy", 15) = {5};
//+
Physical Surface("+yz", 16) = {2};
//+
Physical Surface("-yz", 17) = {1};
//+
Physical Surface("+zx", 18) = {4};
//+
Physical Surface("-zx", 19) = {3};
