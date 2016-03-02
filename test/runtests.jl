using EllipticFEM
using Base.Test

meshpath="mesh_p05.msh";
A = 1.0;
f = (x,y)->5/4*pi^2*cos(x*pi/2)*sin(y*pi);
bddata = [1 2 3 4;
		'D' 'N' 'D' 'N';
		(x,y)->0 (x,y)->-pi/2*sin(x*pi/2)*sin(y*pi) (x,y)->0 (x,y)->pi/2*sin(x*pi/2)*sin(y	*pi)];
mesh,u = solve(meshpath,A,f,bddata);

plot_solution(mesh,u)
