function lphi=lop2d(phi,k,solver)
% function lphi=lop2d(phi,k,solver)

%
% $Id: lop2d.m,v 1.3 2011/03/26 12:56:39 patrick Exp $
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

xa=solver.xa;
xb=solver.xb;
yc=solver.yc;
yd=solver.yd;

nxa=solver.nxa;
nxb=solver.nxb;
nyc=solver.nyc;
nyd=solver.nyd;

nx=solver.nk{k}(1);
ny=solver.nk{k}(2);

xoff=solver.xoff; 
yoff=solver.yoff;

cofx=solver.cofx{k};
cofy=solver.cofy{k};

% Calculate Te on the grid
te=Te2d(solver,k);
% Calculate Ne on the grid
ne=Ne2d(solver,k);

I=1:nx;
J=1:ny;
lphi=phi(I+xoff,J+yoff);

ist=1;
if nxa==1,
	ist=2; 
end;
ifn=nx;
if nxb==1,
	ifn=nx-1;
end;
jst=1;
if nyc==1,
	jst=2; 
end;
jfn=ny;
if nyd==1,
	jfn=ny-1; 
end;

I=ist:ifn;
J=jst:jfn;

Io=I+xoff; Jo=J+yoff;

% L(phi)=d2(phi)/dx2 + d2(phi)/dy2 -exp(phi/Te)
[cofx1,cofy1]=ndgrid(cofx(I,1),cofy(J,1));
[cofx2,cofy2]=ndgrid(cofx(I,2),cofy(J,2));
[cofx3,cofy3]=ndgrid(cofx(I,3),cofy(J,3));

lphi(I,J)= ...
	cofx1.*phi(Io-1,Jo)+cofx2.*phi(Io+1,Jo)+ ...
	cofy1.*phi(Io,Jo-1)+cofy2.*phi(Io,Jo+1)+ ...
	(cofx3+cofy3).*phi(Io,Jo)+nlF(phi(Io,Jo),ne(I,J),te(I,J));
