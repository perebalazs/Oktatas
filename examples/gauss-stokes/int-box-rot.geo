//+
SetFactory("OpenCASCADE");
//+
Box(1) = {0, 0, 0, 2, 2, 2};
//+
Rotate {{0, 1, 0}, {0, 0, 0}, Pi/6} {
  Volume{1}; Point{3}; Point{1}; Point{7}; Point{5}; Point{4}; Point{2}; Point{8}; Point{6}; Curve{2}; Curve{12}; Curve{10}; Curve{6}; Curve{3}; Curve{4}; Curve{1}; Curve{11}; Curve{7}; Curve{9}; Curve{8}; Curve{5}; Surface{6}; Surface{1}; Surface{4}; Surface{3}; Surface{2}; Surface{5}; 
}
//+
Rotate {{1, 0, 0}, {0, 0, 0}, Pi/6} {
  Point{3}; Point{1}; Point{4}; Point{7}; Point{2}; Point{5}; Point{8}; Point{6}; Curve{2}; Curve{3}; Curve{12}; Curve{4}; Curve{1}; Curve{6}; Curve{10}; Curve{7}; Curve{11}; Curve{5}; Curve{8}; Curve{9}; Surface{1}; Surface{6}; Surface{4}; Surface{3}; Surface{2}; Surface{5}; Volume{1}; 
}
//+
Mesh.ElementOrder=2;
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
