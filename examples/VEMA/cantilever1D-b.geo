//+
Rectangle(1) = {0, 0, 0, l, A, 0};
//+
Physical Curve("left", 5) = {4};
//+
Physical Curve("right", 6) = {2};
//+
Physical Curve("bottom", 7) = {1};
//+
Physical Curve("top", 8) = {3};
//+
Physical Surface("body", 9) = {1};

//+
Recombine Surface {1};
Transfinite Line {1,3} = n+1;  // csomópontok száma (elemek száma + 1)
Transfinite Line {2,4} = 2;
Transfinite Surface{1};

Mesh.ElementOrder = order;


Mesh 2;

//+
Point(5) = {0, 0.001, 0, 1.0};
//+
Point(6) = {100, 0.001, 0, 1.0};
//+
Line(5) = {5, 6};
//+
Physical Curve("path", 9) = {5};

