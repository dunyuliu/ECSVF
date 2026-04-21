lcell=1000; 	 // mesh size

pi=3.141592653589793238462643383279502884197169399375;
d1=20000;   	 //depth
dip1=50;    	 //dip_angle
azi1=54.576265;  //azimuth
l1=755306.9;     //fault_length

nline1w=l1/lcell+1;
nline1h=(d1/Cos(pi/180*(90-dip1)))/lcell+1;
origin_x=1096517.000000;
origin_y=4997154.000000;
origin_z=0;
x4=-origin_y-l1*Cos(pi/180*azi1);
y4=origin_x+l1*Sin(pi/180*azi1);
Point(1) = {origin_x,origin_y,origin_z,lcell}; //{0, 0, 0, lcell};
Point(4) = {y4, -x4, origin_z, lcell};

projected_l=d1*Tan(pi/180*(90-dip1));
x2=-origin_y+projected_l*Sin(pi/180*azi1);
y2=origin_x+projected_l*Cos(pi/180*azi1);
x3=x4+projected_l*Sin(pi/180*azi1);
y3=y4+projected_l*Cos(pi/180*azi1);

Point(2) = {y2, -x2, origin_z-d1, lcell};
Point(3) = {y3, -x3, origin_z-d1, lcell};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(1) = {4, 1, 2, 3};
Plane Surface(1) = {1}; //a list of line loop, here loop (1)

Physical Line("wall") = {1,2,3,4};
Transfinite Surface {1} = {4,1,2,3};
Transfinite Line {1,3} = nline1h Using Progression 1;
Transfinite Line {2,4} = nline1w Using Progression 1;
Recombine Surface {1};
Physical Surface("My surface1") = {1};