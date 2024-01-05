//+
SetFactory("OpenCASCADE");
Rectangle(1) = {0, 0, 0, 100, 100, 0};
//+
Physical Curve("left", 5) = {4};
//+
Physical Curve("bottom", 6) = {1};
//+
Physical Curve("top", 7) = {3};
//+
Physical Surface("rect", 8) = {1};
//+
MeshSize {4, 3, 1, 2} = 10;
