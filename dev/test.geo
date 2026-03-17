//+ LowLevelFEM
//+
SetFactory("OpenCASCADE");
//+
Ellipse(1) = {-0.2, -0.3, 0, 0.5, 0.25, 0, 2*Pi};
//+
Curve Loop(1) = {1};
//+
Surface(1) = {1};
