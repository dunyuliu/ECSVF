lcell = 1000;  // mesh size

Point(1) = {1721852.000000, 5421559.000000, -20000.000000, lcell};
Point(2) = {1106147.000000, 4983614.000000, -20000.000000, lcell};
Point(3) = {1096517.000000, 4997154.000000, 0.000000, lcell};
Point(4) = {1712028.000000, 5434959.000000, 0.000000, lcell};

Line(1) = {2, 3};
Line(2) = {3, 4};
Line(3) = {4, 1};
Line(4) = {1, 2};

Line Loop(1) = {1, 2, 3, 4};
Surface(1) = {1};

// Physical groups
Physical Line("wall") = {2, 3, 4, 1};
MeshAlgorithm Surface {1} = {2, 3, 4, 1};
Physical Surface("MySurface1") = {1};