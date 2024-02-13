function [u,f,solver]=fun1(solver)
% function [u,f,solver]=fun1(solver)

%
% $Id: fun1_2d.m,v 1.2 2011/03/26 12:56:38 patrick Exp $
%
% Copyright (c) 2000-2011 Patrick Guio <patrick.guio@gmail.com>
% All Rights Reserved.
%
% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2.  of the License, or (at your
% option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
% Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <http://www.gnu.org/licenses/>.

u=zeros(solver.nx,solver.ny);

global a kx ky X Y phit

a=.01;
kx=10*pi;
ky=8*pi;

fid=fopen('partest1_2d.dat','w');
fprintf(fid,'$\\rho_0=%.2f$, $k_x=%.0f\\pi$, $k_y=%.0f\\pi$ and $T_e=%.0f$',...
	a, kx/pi, ky/pi, solver.temin);
fclose(fid);

x=linspace(0,1,solver.nx);
y=linspace(0,1,solver.ny);
[X,Y]=ndgrid(x,y);

f=ones(solver.nx,solver.ny)+a*cos(kx*X+ky*Y);
phit=cos(kx*X+ky*Y);

