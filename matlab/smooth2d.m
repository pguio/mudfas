function phi=smooth2d(phi,rhs,k,solver)
% function phi=smooth2d(phi,rhs,k,solver)
%
% gauss-seidel-newton red/black point relaxation
%
% L(phi)=d2(phi)/dx2 + d2(phi)/dy2 -exp(phi/Te)
% phi=phi - ( L(phi)-rhs ) / d(L(phi))/dphi

%
% $Id: smooth2d.m,v 1.3 2011/03/26 12:56:39 patrick Exp $
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

% set periodic b.c. indicator
nper=nxa*nxb*nyc*nyd;

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

if nper == 0, % set periodic boundaries if necessary
  phi=per2vb(phi,k,solver);
end

% Calculate Te on the grid
te=Te2d(solver,k);
% Calculate Ne on the grid
ne=Ne2d(solver,k);

[cofx1,cofy1]=ndgrid(cofx(:,1),cofy(:,1));
[cofx2,cofy2]=ndgrid(cofx(:,2),cofy(:,2));
[cofx3,cofy3]=ndgrid(cofx(:,3),cofy(:,3));

% L(phi)=d2(phi)/dx2 + d2(phi)/dy2 -exp(phi/Te)
% phi=phi - ( L(phi)-rhs ) / d(L(phi))/dphi

% relax on red grid points
I=ist:2:ifn;
J=jst:2:jfn;
Io=I+xoff; Jo=J+yoff;
phi=Gauss_Seidel_Newton_2d(phi,rhs, ...
	cofx1,cofx2,cofx3,cofy1,cofy2,cofy3,ne,te,I,J,Io,Jo);

I=2:2:ifn;
J=2:2:jfn;
Io=I+xoff; Jo=J+yoff;
phi=Gauss_Seidel_Newton_2d(phi,rhs, ...
	cofx1,cofx2,cofx3,cofy1,cofy2,cofy3,ne,te,I,J,Io,Jo);

if nper == 0, 
	  phi=per2vb(phi,k,solver);
end

% relax on black grid points
I=ist:2:ifn;
J=2:2:jfn;
Io=I+xoff; Jo=J+yoff;
phi=Gauss_Seidel_Newton_2d(phi,rhs, ...
	cofx1,cofx2,cofx3,cofy1,cofy2,cofy3,ne,te,I,J,Io,Jo);

I=2:2:ifn;
J=jst:2:jfn;
Io=I+xoff; Jo=J+yoff;
phi=Gauss_Seidel_Newton_2d(phi,rhs, ...
	cofx1,cofx2,cofx3,cofy1,cofy2,cofy3,ne,te,I,J,Io,Jo);
	
if nper == 0, % set periodic boundaries
	phi=per2vb(phi,k,solver);
end

function phi=per2vb(phi,k,solver)
% set virtual periodic boundaries from interior values
% in two dimensions (for all 2-d solvers)

nxa=solver.nxa;
nxb=solver.nxb;
nyc=solver.nyc;
nyd=solver.nyd;

nx=solver.nk{k}(1);
ny=solver.nk{k}(2);

xoff=solver.xoff;
yoff=solver.yoff;

if nxa*nxb==0,
	J=1:ny;
	phi(0+xoff,J+yoff)=phi(nx-1+xoff,J+yoff);
	phi(nx+xoff,J+yoff)=phi(1+xoff,J+yoff);
	phi(nx+1+xoff,J+yoff)=phi(2+xoff,J+yoff);
end

if nyc*nyd==0,
	I=1:nx;
	phi(I+xoff,0+yoff)=phi(I+xoff,ny-1+yoff);
	phi(I+xoff,ny+yoff)=phi(I+xoff,1+yoff);
	phi(I+xoff,ny+1+yoff)=phi(I+xoff,2+yoff);
end

function phi=Gauss_Seidel_Newton_2d(phi,rhs, ...
	cofx1,cofx2,cofx3,cofy1,cofy2,cofy3,ne,te,I,J,Io,Jo)
%

phi(Io,Jo)=phi(Io,Jo)-( ...
	cofx1(I,J).*phi(Io-1,Jo)+cofx2(I,J).*phi(Io+1,Jo)+ ...
	cofy1(I,J).*phi(Io,Jo-1)+cofy2(I,J).*phi(Io,Jo+1)+ ...
	(cofx3(I,J)+cofy3(I,J)).*phi(Io,Jo)+ ...
	nlF(phi(Io,Jo),ne(I,J),te(I,J))-rhs(I,J))./ ...
	(cofx3(I,J)+cofy3(I,J)+nldF(phi(Io,Jo),ne(I,J),te(I,J)));
