function plotEfield3d(phi,solver,argcaxis,mnp,func)
% function plotEfield3d(phi,solver,argcaxis,mnp,func)

%
% $Id: plotEfield3d.m,v 1.2 2011/03/26 12:56:38 patrick Exp $
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

hx=(solver.xb-solver.xa)/(solver.nx-1);
hy=(solver.yd-solver.yc)/(solver.ny-1);
hz=(solver.zf-solver.ze)/(solver.nz-1);

m=mnp(1); n=mnp(2); p=mnp(3);

[Ex,Ey,Ez]=gradient(phi,hx,hy,hz);

subplot(m,n,p), imagesc(x,y,Ex(:,:,fix(end/2))); colorbar('v');
if mod(p-1,n)~=0, 
  set(gca,'yticklabel',[]); 
end;
if fix((p-1)/n)~=m-1, 
  set(gca,'xticklabel',[]); 
end;
title(['E_x ' func]);

subplot(m,n,p+1), imagesc(x,y,Ey(:,:,fix(end/2))); 
if mod(p,n)~=0, set(gca,'yticklabel',[]); end;
if fix(p/n)~=m-1, set(gca,'xticklabel',[]); end;
caxis(argcaxis); colorbar('v');
title(['E_y ' func]);

subplot(m,n,p+2), imagesc(x,y,Ez(:,:,fix(end/2))); 
if mod(p,n)~=0, set(gca,'yticklabel',[]); end;
if fix(p/n)~=m-1, set(gca,'xticklabel',[]); end;
caxis(argcaxis); colorbar('v');
title(['E_y ' func]);

drawnow

