function plotphi3d(phi,solver,argcaxis,mnp,func)
% function plotphi3d(phi,solver,argcaxis,mnp,func)

%
% $Id: plotphi3d.m,v 1.2 2011/03/26 12:56:38 patrick Exp $
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

m=mnp(1); n=mnp(2); p=mnp(3);

subplot(m,n,p), 
imagesc(x,y,phi(:,:,fix(end/2))'); 
xlabel('x')
ylabel('y')
axis xy
colorbar('v');
if ~isempty(argcaxis),
	caxis(argcaxis); colorbar
end
if mod(p-1,n)~=0, set(gca,'yticklabel',[]); end;
%if fix((p-1)/n)~=m-1, set(gca,'xticklabel',[]); end;
title(['\phi' func]);

drawnow
