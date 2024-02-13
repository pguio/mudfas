function [u,f,solver]=fun6(solver)
% function [u,f,solver]=fun6(solver)

%
% $Id: fun6_3d.m,v 1.2 2011/03/26 12:56:39 patrick Exp $
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


u=zeros(solver.nx,solver.ny, solver.nz);


randn('state',0);
f=ones(solver.nx,solver.ny,solver.nz)+0.0*randn(solver.nx,solver.ny,solver.nz);
a=.05;
kx=10*pi;
ky=8*pi;
ky=0*pi;


solver_const;

solver.temax=solver.temin;
solver.teshape=Te3dModel.Const3d;

solver.nemin=1.0;
solver.nemax=1.1;
solver.neshape=Ne3dModel.XY_Gauss;

[X,Y,Z]=ndgrid(linspace(-1,1,solver.nx),linspace(-1,1,solver.ny),zeros(1,solver.nz));
f=(solver.nemax-solver.nemin)*exp(-(X.^2+Y.^2)/2/solver.sigman^2)+solver.nemin;

x=linspace(0,1,solver.nx);
y=linspace(0,1,solver.ny);
z=linspace(0,1,solver.nz);
[X,Y,Z]=ndgrid(x,y,z);
f=f+a*cos(kx*X+ky*Y)+a/2*randn(solver.nx,solver.ny,solver.nz);
