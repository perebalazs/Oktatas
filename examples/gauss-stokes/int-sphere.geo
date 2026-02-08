//+
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 10, 0, 1.0};
//+
Point(2) = {0, 0, 10, 1.0};
//+
Point(3) = {0, 0, 0, 1.0};
//+
Circle(1) = {1, 3, 2};
//+
Extrude {{0, 1, 0}, {0, 0, 0}, Pi/2} {
  Curve{1}; 
}
MeshSize {:} = 10;
Mesh.ElementOrder=2;
Mesh 2;
//+
Physical Surface("surf", 5) = {1};
//+
Physical Curve("x", 6) = {1};
//+
Physical Curve("y", 7) = {3};
//+
Physical Curve("z", 8) = {4};
