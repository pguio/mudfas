function plotdefect2d(u,f,solver)
% function plotdefect2d(u,f,solver)

%
% $Id: plotdefect2d.m,v 1.4 2011/03/26 12:56:39 patrick Exp $
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

x=linspace(solver.xa,solver.xb,solver.nx);
y=linspace(solver.yc,solver.yd,solver.ny);

subplot(121)
imagesc(x,y,u')
axis xy
title('\phi')
colorbar('h')

subplot(122)

%d=lop2d(u{kgrid},f{kgrid},kgrid,solver)-f{kgrid};
hx=(solver.xb-solver.xa)/(solver.nx-1);
hy=(solver.yd-solver.yc)/(solver.ny-1);

U=zeros(solver.nx+2,solver.ny+2);
U(2:end-1,2:end-1)=u;
U([1 end],2:end-1)=u([1 end],:);
U(2:end-1,[1 end])=u(:,[1 end]);
U([1 end],[1 end])=u([1 end],[1 end]);
L=4.0*del2(U,hy,hx);

d=L(2:end-1,2:end-1)+ ...
	nlF(u,Ne2d(solver,solver.ngrid),Te2d(solver,solver.ngrid))-f;
d=d(2:end-1,2:end-1); 
x=x(2:end-1); 
y=y(2:end-1);
imagesc(x,y,d')
title('defect')
axis xy,
colorbar('h'),


drawnow

fprintf(1,'mean |defect|=%.4e std |defect|=%.4e\n', ...
	mean(abs(d(:))), std(abs(d(:))));

