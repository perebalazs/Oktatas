//+
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
l = DefineNumber[ 100, Name "Parameters/l" ];
//+
n = DefineNumber[ 2, Name "Parameters/n" ];
//+
F = DefineNumber[ 10, Name "Parameters/F" ];
//+
fx = DefineNumber[ 1, Name "Parameters/fx" ];
//+
Point(2) = {l, 0, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Transfinite Curve {1} = n+1 Using Progression 1;
//+
Mesh 1;
//+
Physical Point("left", 2) = {1};
//+
Physical Point("right", 3) = {2};
//+
Physical Curve("truss", 4) = {1};
