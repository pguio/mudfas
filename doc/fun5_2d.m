function [u,f,solver]=fun5(solver)
% function [u,f,solver]=fun5(solver)

%
% $Id: fun5_2d.m,v 1.2 2011/03/26 12:56:38 patrick Exp $
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

global nx ny

solver.nx=nx;
solver.ny=ny;

solver.xa=0;
solver.xb=solver.nx-1;
solver.yc=0;
solver.yd=solver.ny-1;

u=zeros(solver.nx,solver.ny);

global a kx ky X Y x y phit

a=.01;
kx=0.0;
ky=12*pi;

fid=fopen('partest5_2d.dat','w');
fprintf(fid,'$k_y=%.0f\\pi$', ky/pi);
fclose(fid);

x=linspace(0,1,solver.nx);
y=linspace(0,1,solver.ny);
[X,Y]=ndgrid(x,y);

f=ones(solver.nx,solver.ny)+a*cos(kx*X+ky*Y);
phit=cos(kx*X+ky*Y);
