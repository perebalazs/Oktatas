//+
SetFactory("OpenCASCADE");
Disk(1) = {0, 15, 0, 5, 5};

Mesh 2;
//+
Physical Surface("disk", 2) = {1};
