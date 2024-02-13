function phi=smooth3d(phi,rhs,k,solver)
% function phi=smooth3d(phi,rhs,k,solver)
%
% gauss-seidel point relaxation with red/black ordering
% in three dimensions for nonseparable pde
% relax in order:
% (1) red (x,y) on odd z planes
% (2) black (x,y) on even z planes
% (3) black (x,y) on odd z planes
% (4) red (x,y) on even z planes
%
% L(phi)=d2(phi)/dx2 + d2(phi)/dy2 + d3(phi)/dy3 -exp(phi/Te)
% phi=phi - ( L(phi)-rhs ) / d(L(phi))/dphi

%
% $Id: smooth3d.m,v 1.3 2011/03/26 12:56:39 patrick Exp $
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

% set periodic b.c. indicator
nper=nxa*nxb*nyc*nyd*nze*nzf;

% set loop limits to avoid specified boundaries
% in red/black sweeps
ist=1;
if nxa==1, 
	ist=3; 
end
ifn=nx;
if nxb==1, 
	ifn=nx-1; 
end
jst=1;
if nyc==1, 
	jst=3; 
end
jfn=ny;
if nyd==1, 
	jfn=ny-1; 
end
kst=1;
if nze==1,
	kst=3;
end;
kfn=nz;
if nzf==1,
	kfn=nz-1;
end

if nper == 0, % set periodic boundaries if necessary
	phi=per3vb(phi,k,solver);
end

% Calculate Te on the grid
te=Te3d(solver,k);
% Calculate Ne on the grid
ne=Ne3d(solver,k);

[cofx1,cofy1,cofz1]=ndgrid(cofx(:,1),cofy(:,1),cofz(:,1));
[cofx2,cofy2,cofz2]=ndgrid(cofx(:,2),cofy(:,2),cofz(:,2));
[cofx3,cofy3,cofz3]=ndgrid(cofx(:,3),cofy(:,3),cofz(:,3));


% red (x,y) on odd z planes
I=ist:2:ifn;
J=jst:2:jfn;
K=kst:2:kfn;
Io=I+xoff; Jo=J+yoff; Ko=K+zoff;
phi=Gauss_Seidel_Newton_2d(phi,rhs, ...
	cofx1,cofx2,cofx3,cofy1,cofy2,cofy3,cofz1,cofz2,cofz3, ...
	ne,te,I,J,K,Io,Jo,Ko);

I=2:2:ifn;
J=2:2:jfn;
Io=I+xoff; Jo=J+yoff;
phi=Gauss_Seidel_Newton_2d(phi,rhs, ...
	cofx1,cofx2,cofx3,cofy1,cofy2,cofy3,cofz1,cofz2,cofz3, ...
	ne,te,I,J,K,Io,Jo,Ko);

if nper == 0, 
	phi=per3vb(phi,k,solver);
end

% black (x,y) or even z planes
I=ist:2:ifn;
K=2:2:kfn;
Io=I+xoff; Ko=K+zoff;
phi=Gauss_Seidel_Newton_2d(phi,rhs, ...
	cofx1,cofx2,cofx3,cofy1,cofy2,cofy3,cofz1,cofz2,cofz3, ...
	ne,te,I,J,K,Io,Jo,Ko);

I=2:2:ifn;
J=jst:2:jfn;
Io=I+xoff; Jo=J+yoff;
phi=Gauss_Seidel_Newton_2d(phi,rhs, ...
	cofx1,cofx2,cofx3,cofy1,cofy2,cofy3,cofz1,cofz2,cofz3, ...
	ne,te,I,J,K,Io,Jo,Ko);

if nper == 0,
	phi=per3vb(phi,k,solver);
end

% black (x,y) on odd z planes
I=ist:2:ifn;
J=2:2:jfn;
K=kst:2:kfn;
Io=I+xoff; Jo=J+yoff; Ko=K+zoff;
phi=Gauss_Seidel_Newton_2d(phi,rhs, ...
	cofx1,cofx2,cofx3,cofy1,cofy2,cofy3,cofz1,cofz2,cofz3, ...
	ne,te,I,J,K,Io,Jo,Ko);

I=2:2:ifn;
J=jst:2:jfn;
Io=I+xoff; Jo=J+yoff;
phi=Gauss_Seidel_Newton_2d(phi,rhs, ...
	cofx1,cofx2,cofx3,cofy1,cofy2,cofy3,cofz1,cofz2,cofz3, ...
	ne,te,I,J,K,Io,Jo,Ko);

if nper == 0,
	phi=per3vb(phi,k,solver);
end

% red(x,y) on even z planes
I=ist:2:ifn;
K=2:2:kfn;
Io=I+xoff; Ko=K+zoff;
phi=Gauss_Seidel_Newton_2d(phi,rhs, ...
	cofx1,cofx2,cofx3,cofy1,cofy2,cofy3,cofz1,cofz2,cofz3, ...
	ne,te,I,J,K,Io,Jo,Ko);

I=2:2:ifn;
J=2:2:jfn;
Io=I+xoff; Jo=J+yoff;
phi=Gauss_Seidel_Newton_2d(phi,rhs, ...
	cofx1,cofx2,cofx3,cofy1,cofy2,cofy3,cofz1,cofz2,cofz3, ...
	ne,te,I,J,K,Io,Jo,Ko);

if nper == 0, % set periodic boundaries if necessary
	phi=per3vb(phi,k,solver);
end

function phi=per3vb(phi,k,solver)
% set virtual periodic boundaries from interior values
% in three dimensions (for all 3-d solvers)

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

if nxa*nxb==0,
	K=1:nz;
	J=1:ny;
	phi(0+xoff,J+yoff,K+zoff) = phi(nx-1+xoff,J+yoff,K+zoff);
	phi(nx+xoff,J+yoff,K+zoff) = phi(1+xoff,J+yoff,K+zoff);
	phi(nx+1+xoff,J+yoff,K+zoff) = phi(2+xoff,J+yoff,K+zoff);
end

if nyc*nyd==0,
	K=1:nz;
	I=1:nx;
	phi(I+xoff,0+yoff,K+zoff) = phi(I+xoff,ny-1+yoff,K+zoff);
	phi(I+xoff,ny+yoff,K+zoff) = phi(I+xoff,1+yoff,K+zoff);
	phi(I+xoff,ny+1+yoff,K+zoff) = phi(I+xoff,2+yoff,K+zoff);
end

if nze*nzf==0,
	J=1:ny;
	I=1:nx;
	phi(I+xoff,J+yoff,0+zoff) = phi(I+xoff,J+yoff,nz-1+zoff);
	phi(I+xoff,J+yoff,nz+zoff) = phi(I+xoff,J+yoff,1+zoff);
	phi(I+xoff,J+yoff,nz+1+zoff) = phi(I+xoff,J+yoff,2+zoff);
end

function phi=Gauss_Seidel_Newton_2d(phi,rhs, ...
	cofx1,cofx2,cofx3,cofy1,cofy2,cofy3,cofz1,cofz2,cofz3, ...
	ne,te,I,J,K,Io,Jo,Ko)

% L(phi)=d2(phi)/dx2 + d2(phi)/dy2 + d3(phi)/dy3 -exp(phi/Te)
% phi=phi - ( L(phi)-rhs ) / d(L(phi))/dphi

phi(Io,Jo,Ko)=phi(Io,Jo,Ko)-( ...
	cofx1(I,J,K).*phi(Io-1,Jo,Ko)+cofx2(I,J,K).*phi(Io+1,Jo,Ko)+ ...
	cofy1(I,J,K).*phi(Io,Jo-1,Ko)+cofy2(I,J,K).*phi(Io,Jo+1,Ko)+ ...
	cofz1(I,J,K).*phi(Io,Jo,Ko-1)+cofz2(I,J,K).*phi(Io,Jo,Ko+1)+ ...
	(cofx3(I,J,K)+cofy3(I,J,K)+cofz3(I,J,K)).*phi(Io,Jo,Ko)+ ...
	nlF(phi(Io,Jo,Ko),ne(I,J,K),te(I,J,K))-rhs(I,J,K))./( ...
	cofx3(I,J,K)+cofy3(I,J,K)+cofz3(I,J,K)+ ...
	nldF(phi(Io,Jo,Ko),ne(I,J,K),te(I,J,K)));

