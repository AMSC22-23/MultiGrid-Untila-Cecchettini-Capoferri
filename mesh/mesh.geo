Lx = 2.0;
Ly = 2.0;


h = 1.0;

Point(1) = {0, 0, 0, h};
Point(2) = {Lx, 0, 0, h};
Point(3) = {Lx, Ly, 0, h};
Point(4) = {0, Ly, 0, h};

Line(5) = {1, 2};
Line(6) = {2, 3};
Line(7) = {3, 4};
Line(8) = {4, 1};

Line Loop(1) = {5, 6, 7, 8};

Plane Surface(1) = {1};

Physical Surface(10) = {1};

Mesh 2;



