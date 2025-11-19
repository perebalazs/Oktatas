//+
SetFactory("OpenCASCADE");
//+
R = 15; // > r
r = 5;
//+
a = 1.7;
//+
Circle(1) = {0, R, 0, r, 0, 2*Pi};
//+
Point(2) = {-a/5*r, R-a/5*r, 0, 1.0};
//+
Point(3) = {a/5*r, R-a/5*r, 0, 1.0};
//+
Point(4) = {a/5*r, R+a/5*r, 0, 1.0};
//+
Point(5) = {-a/5*r, R+a/5*r, 0, 1.0};
//+
Point(6) = {-r, R-r, 0, 1.0};
//+
Point(7) = {r, R-r, 0, 1.0};
//+
Point(8) = {r, R+r, 0, 1.0};
//+
Point(9) = {-r, R+r, 0, 1.0};
//+
Line(2) = {5, 9};
//+
Line(3) = {4, 8};
//+
Line(4) = {3, 7};
//+
Line(5) = {2, 6};
//+
Line(6) = {2, 5};
//+
Line(7) = {5, 4};
//+
Line(8) = {4, 3};
//+
Line(9) = {3, 2};
//+
Curve Loop(1) = {1};
//+
Plane Surface(1) = {1};
//+
BooleanFragments{ Surface{1}; Delete; }{ Curve{3}; Curve{8}; Curve{7}; Curve{4}; Curve{9}; Curve{6}; Curve{2}; Curve{5}; Delete; }
//+
Recursive Delete {
  Curve{21}; Curve{19}; Curve{20}; Curve{22}; 
}
//+
Coherence;
//+
pcs1 = 6;
pcs2 = Ceil(((R-r)*Pi/3) / (2*r*Pi / (pcs1*8)));
Transfinite Curve {10, 13} = pcs1+1 Using Progression 1;
//+
Transfinite Curve {16, 18, 15, 7, 8, 9, 6} = pcs1*2+1 Using Progression 1;
//+
Transfinite Curve {17, 11, 12, 14} = pcs1*2+1 Using Progression 1;
//+
Transfinite Surface {4} = {10, 2, 5, 9};
//+
Transfinite Surface {2} = {5, 4, 8, 9};
//+
Transfinite Surface {1} = {3, 6, 8, 4};
//+
Transfinite Surface {3} = {10, 6, 3, 2};
//+
Transfinite Surface {5} = {2, 3, 4, 5};
//+
Recombine Surface {4, 5, 2, 3, 1};

Extrude {{1, 0, 0}, {0, 0, 0}, Pi/4} {
  Surface{3}; Surface{1}; Surface{2}; Surface{4}; Surface{5}; Layers{pcs2+1}; Recombine;
}
//+
Extrude {{1, 0, 0}, {0, 0, 0}, Pi/4} {
  Surface{23}; Surface{22}; Surface{19}; Surface{15}; Surface{10}; Layers{pcs2+1}; Recombine;
}

Coherence;
Mesh.ElementOrder=1;
Mesh 3;
//+
Physical Volume("body", 63) = {4, 3, 5, 1, 7, 2, 8, 10, 6, 9};
//+
Physical Surface("supp", 64) = {2, 1, 3, 4, 5};
//+
Physical Surface("load", 65) = {41, 32, 35, 39, 28};
//+
Point(29) = {0, (R+r)*Sqrt(2)/2, (R+r)*Sqrt(2)/2, 1.0};
//+
Point(30) = {0, (R-r)*Sqrt(2)/2, (R-r)*Sqrt(2)/2, 1.0};
//+
Line(63) = {29, 30};
//+
Physical Curve("path", 66) = {63};
//+
Physical Surface("mid", 67) = {22, 10, 15, 19, 23};
//+
Physical Surface("norm1", 68) = {34};
