//+
SetFactory("OpenCASCADE");
//+
Point(1) = {-0.7, 0.4, 0, 1.0};
//+
Point(2) = {0.1, -0, 0, 1.0};
//+
Point(3) = {-0.4, -0.4, 0, 1.0};
//+
Point(4) = {-0.8, -0.2, 0, 1.0};
//+
Line(1) = {4, 3};
//+
Line(2) = {2, 3};
//+
Line(3) = {2, 1};
//+
Line(4) = {1, 4};
//+
Extrude {0, 0, 1} {
  Curve{4}; Curve{3}; Curve{2}; Curve{1}; 
}
//+
Curve Loop(5) = {7, 12, 11, 9};
//+
Surface(5) = {5};
