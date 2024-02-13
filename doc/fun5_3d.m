function [u,f,solver]=fun5(solver)
% function [u,f,solver]=fun5(solver)

%
% $Id: fun5_3d.m,v 1.2 2011/03/26 12:56:38 patrick Exp $
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

global nx ny nz

solver.nx=nx;
solver.ny=ny;
solver.nz=nz;

solver.xa=0;
solver.xb=solver.nx-1;
solver.yc=0;
solver.yd=solver.ny-1;
solver.ze=0;
solver.zf=solver.nz-1;

u=zeros(solver.nx,solver.ny,solver.nz);

global a kx ky kz X Y Z x y z phit

a=.01;
kx=0.0;
ky=12*pi;
kz=0.0;

fid=fopen('partest5_3d.dat','w');
fprintf(fid,'$k_y=%.0f\\pi$', ky/pi);
fclose(fid);

x=linspace(0,1,solver.nx);
y=linspace(0,1,solver.ny);
z=linspace(0,1,solver.nz);
[X,Y,Z]=ndgrid(x,y,z);

f=ones(solver.nx,solver.ny,solver.nz)+a*cos(kx*X+ky*Y+kz*Y);
phit=cos(kx*X+ky*Y+kz*Y);
