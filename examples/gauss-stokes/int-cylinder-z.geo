//+
SetFactory("OpenCASCADE");
//+
Point(1) = {10, 0, 0, 1.0};
//+
Point(2) = {0, 10, 0, 1.0};
//+
Point(3) = {0, 0, 0, 1.0};
//+
Circle(1) = {1, 3, 2};
//+
Extrude {0, 0, 10} {
  Curve{1}; 
}

Mesh 2;
//+
Physical Surface("surf", 5) = {1};
