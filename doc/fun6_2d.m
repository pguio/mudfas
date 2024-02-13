function [u,f,solver]=fun6(solver)
% function [u,f,solver]=fun6(solver)

%
% $Id: fun6_2d.m,v 1.2 2011/03/26 12:56:38 patrick Exp $
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

global a kx ky X Y phit

x=linspace(0,1,solver.nx);
y=linspace(0,1,solver.ny);
[X,Y]=ndgrid(x,y);

nb=1.2;
f=ones(solver.nx,solver.ny);
imn=fix(solver.nx/3);
imx=solver.nx-imn;
jmn=fix(solver.ny/2);
jmx=solver.ny;
f(imn:imx,jmn:jmx)=nb*f(imn:imx,jmn:jmx);

randn('state',0);
m=0.0;
sd=0.0;
f=f+sd*(randn(solver.nx,solver.ny)-m);

phit=zeros(size(f));
phit([imn:imx],[jmn:jmx])=solver.temin*log(nb);

fid=fopen('partest6_2d.dat','w');
fprintf(fid,'$n_b=%.1f$ and $T_e=%.0f$',nb, solver.temin);
fclose(fid);

