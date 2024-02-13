function [u,f,solver]=fun3(solver)
% function [u,f,solver]=fun3(solver)

%
% $Id: fun3_2d.m,v 1.2 2011/03/26 12:56:39 patrick Exp $
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

f=ones(solver.nx,solver.ny);

a=0.05;
kx=10*pi;
ky=8*pi;

x=linspace(0,1,solver.nx);
y=linspace(0,1,solver.ny);
[X,Y]=ndgrid(x,y);

f=f+a*cos(kx*X+ky*Y);

imn=fix(solver.nx/3); 
imx=solver.nx;
jmn=fix(solver.ny/3); 
jmx=fix(2*solver.ny/3);

f(imn:imx,jmn:jmx)=1.2*f(imn:imx,jmn:jmx);


