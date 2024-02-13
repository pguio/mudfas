function lphi=lop3d(phi,k,solver)
% function lphi=lop3d(phi,k,solver)

%
% $Id: lop3d.m,v 1.3 2011/03/26 12:56:39 patrick Exp $
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
ze=solver.ze;
zf=solver.zf;

nxa=solver.nxa;
nxb=solver.nxb;
nyc=solver.nyc;
nyd=solver.nyd;
nze=solver.nze;
nzf=solver.nzf;

nx=solver.nk{k}(1);
ny=solver.nk{k}(2);
nz=solver.nk{k}(3);

xoff=solver.xoff;
yoff=solver.yoff;
zoff=solver.zoff;

cofx=solver.cofx{k};
cofy=solver.cofy{k};
cofz=solver.cofz{k};

% calculate Te on the grid k
te=Te3d(solver,k);
% calculate Ne on the grid k
ne=Ne3d(solver,k);

I=1:nx;
J=1:ny;
K=1:nz;
lphi=phi(I+xoff,J+yoff,K+zoff);

% set loop limits
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
kst=1;
if nze==1,
	kst=2; 
end;
kfn=nz;
if nzf==1,
	kfn=nz-1; 
end;

% compute fine grid residual
I=ist:ifn;
J=jst:jfn;
K=kst:kfn;

Io=I+xoff; 
Jo=J+yoff; 
Ko=K+zoff;

[cofx1,cofy1,cofz1]=ndgrid(cofx(I,1),cofy(J,1),cofz(K,1));
[cofx2,cofy2,cofz2]=ndgrid(cofx(I,2),cofy(J,2),cofz(K,2));
[cofx3,cofy3,cofz3]=ndgrid(cofx(I,3),cofy(J,3),cofz(K,3));

% L(phi)=d2(phi)/dx2 + d2(phi)/dy2 + d2(phi)/dz2 -exp(phi/Te)
lphi(I,J,K)= ...
	cofx1.*phi(Io-1,Jo,Ko)+cofx2.*phi(Io+1,Jo,Ko)+ ...
	cofy1.*phi(Io,Jo-1,Ko)+cofy2.*phi(Io,Jo+1,Ko)+ ...
	cofz1.*phi(Io,Jo,Ko-1)+cofz2.*phi(Io,Jo,Ko+1)+ ...
	(cofx3+cofy3+cofz3).*phi(Io,Jo,Ko)+nlF(phi(Io,Jo,Ko),ne(I,J,K),te(I,J,K));
