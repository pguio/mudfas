function plotdefect3d(u,f,solver)
% function plotdefect3d(u,f,solver)

%
% $Id: plotdefect3d.m,v 1.4 2011/03/26 12:56:39 patrick Exp $
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
z=linspace(solver.ze,solver.zf,solver.nz);

subplot(321)
imagesc(x,y,squeeze(u(:,:,fix(end/2)))')
xlabel('x')
ylabel('y')
title('\phi')
axis xy
colorbar('v')

subplot(323)
imagesc(x,z,squeeze(u(:,fix(end/2),:))')
xlabel('x')
ylabel('z')
title('\phi')
axis xy
colorbar('v')

subplot(325)
imagesc(y,z,squeeze(u(fix(end/2),:,:))')
xlabel('y')
ylabel('z')
title('\phi')
axis xy
colorbar('v')

% d=lop3d(u{K},f{K},K,solver)-f{K};
hx=(solver.xb-solver.xa)/(solver.nx-1);
hy=(solver.yd-solver.yc)/(solver.ny-1);
hz=(solver.zf-solver.ze)/(solver.nz-1);

U=zeros(solver.nx+2,solver.ny+2,solver.nz+2);
U(2:end-1,2:end-1,2:end-1)=u;
U([1 end],2:end-1,2:end-1)=u([1 end],:,:);
U(2:end-1,[1 end],2:end-1)=u(:,[1 end],:);
U(2:end-1,2:end-1,[1 end])=u(:,:,[1 end]);
U([1 end],[1 end],[1 end])=u([1 end],[1 end],[1 end]);
L=6.0*del2(U,hy,hx,hz);

d=L(2:end-1,2:end-1,2:end-1)+ ...
	nlF(u,Ne3d(solver,solver.ngrid),Te3d(solver,solver.ngrid))-f;
d=d(2:end-1,2:end-1,2:end-1); 
x=x(2:end-1); 
y=y(2:end-1); 
z=z(2:end-1);

subplot(322)
imagesc(x,y,squeeze(d(:,:,fix(end/2)))')
xlabel('x')
ylabel('y')
title('defect')
axis xy, 
colorbar('v'),

subplot(324)
imagesc(x,z,squeeze(d(:,fix(end/2),:))')
xlabel('x')
ylabel('z')
title('defect')
axis xy, 
colorbar('v'),

subplot(326)
imagesc(y,z,squeeze(d(fix(end/2),:,:))')
xlabel('y')
ylabel('z')
title('defect')
axis xy, 
colorbar('v'),

drawnow

fprintf(1,'mean |defect|=%.4e std |defect|=%.4e\n', ...
  mean(abs(d(:))), std(abs(d(:))));

