function rhsc=res3(resf,rhsc,k,solver)
% function rhsc=res3(resf,rhsc,k,solver)
%
% restrict resf grid residual in resf to coarse grid in rhsc
% using full weighting for all double precision 3d codes

%
% $Id: res3.m,v 1.3 2011/03/26 12:56:39 patrick Exp $
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

nxa=solver.nxa;
nxb=solver.nxb;
nyc=solver.nyc;
nyd=solver.nyd;
nze=solver.nze;
nzf=solver.nzf;

nx=solver.nk{k+1}(1);
ny=solver.nk{k+1}(2);
nz=solver.nk{k+1}(3);

ncx=solver.nk{k}(1);
ncy=solver.nk{k}(2);
ncz=solver.nk{k}(3);

% set x,y coarsening integer subscript scales
ix=1;
if ncx==nx,
	ix=0; 
end;
jy=1;
if ncy==ny,
	jy=0;
end;
kz=1;
if ncz==nz,
	kz=0;
end;

% restrict on interior

IC=2:ncx-1;
I=IC+ix*(IC-1);
JC=2:ncy-1;
J=JC+jy*(JC-1);
KC=2:ncz-1;
K=KC+kz*(KC-1);

IM1=I-1; IP1=I+1;
JM1=J-1; JP1=J+1;
KM1=K-1; KP1=K+1;
rhsc=Full_Weighting_3d(resf,IM1,I,IP1,JM1,J,JP1,KM1,K,KP1,rhsc,IC,JC,KC);

% set residual on boundaries

% x=xa and x=xb
IC=[1 ncx];
I=IC+ix*(IC-1);
IM1=max(I-1,2);
IP1=min(I+1,nx-1);
IM1(find(I==1 & nxa==0))=nx-1;
IP1(find(I==nx & nxb==0))=2;

% (y,z) interior
JC=2:ncy-1;
J=JC+jy*(JC-1);
KC=2:ncz-1;
K=KC+kz*(KC-1);

JM1=J-1; JP1=J+1;
KM1=K-1; KP1=K+1;
rhsc=Full_Weighting_3d(resf,IM1,I,IP1,JM1,J,JP1,KM1,K,KP1,rhsc,IC,JC,KC);

% x=xa,xb and y=yc,yd interior edges
JC=[1 ncy];
J=JC+jy*(JC-1);
JM1=max(J-1,2);
JP1=min(J+1,ny-1);
JM1(find(J==1 & nyc==0))=ny-1;
JP1(find(J==ny & nyd==0))=2;

rhsc=Full_Weighting_3d(resf,IM1,I,IP1,JM1,J,JP1,KM1,K,KP1,rhsc,IC,JC,KC);

% x=xa,xb; y=yc,yd; z=ze,zf corners
KC=[1 ncz];
K=KC+kz*(KC-1);
KM1=max(K-1,2);
KP1=min(K+1,nz-1);
KM1(find(K==1 & nze==0))=nz-1;
KP1(find(K==nz & nzf==0))=2;

rhsc=Full_Weighting_3d(resf,IM1,I,IP1,JM1,J,JP1,KM1,K,KP1,rhsc,IC,JC,KC);

% x=xa,xb and z=ze,zf edges
JC=2:ncy-1;
J=JC+jy*(JC-1);

JM1=J-1; JP1=J+1;
rhsc=Full_Weighting_3d(resf,IM1,I,IP1,JM1,J,JP1,KM1,K,KP1,rhsc,IC,JC,KC);

% y boundaries y=yc and y=yd
JC=[1 ncy];
J=JC+jy*(JC-1);
JM1=max(J-1,2);
JP1=min(J+1,ny-1);
JM1(find(J==1 & nyc==0))=ny-1;
JP1(find(J==ny & nyd==0))=2;

% (x,z) interior
IC=2:ncx-1;
I=IC+ix*(IC-1);
KC=2:ncz-1;
K=KC+kz*(KC-1);

IM1=I-1; IP1=I+1;
KM1=K-1; KP1=K+1;
rhsc=Full_Weighting_3d(resf,IM1,I,IP1,JM1,J,JP1,KM1,K,KP1,rhsc,IC,JC,KC);

% y=yc,yd and z=ze,zf edges
KC=[1 ncz];
K=KC+kz*(KC-1);
KM1=max(K-1,2);
KP1=min(K+1,nz-1);
KM1(find(K==1 & nze==0))=nz-1;
KP1(find(K==nz & nzf==0))=2;

% interior in x
rhsc=Full_Weighting_3d(resf,IM1,I,IP1,JM1,J,JP1,KM1,K,KP1,rhsc,IC,JC,KC);

% (x,y) interior
JC=2:ncy-1;
J=JC+jy*(JC-1);

JM1=J-1; JP1=J+1;
rhsc=Full_Weighting_3d(resf,IM1,I,IP1,JM1,J,JP1,KM1,K,KP1,rhsc,IC,JC,KC);

return;

% set coarse grid residual to zero at specified boundaries
if nxa==1,
	IC=1;
	KC=1:ncz;
	JC=1:ncy;
	rhsc(IC,JC,KC) = 0.0;
end
if nxb==1,
	IC=ncx;
	KC=1:ncz;
	JC=1:ncy;
	rhsc(IC,JC,KC) = 0.0;
end
if nyc==1,
	JC=1;
	KC=1:ncz;
	IC=1:ncx;
	rhsc(IC,JC,KC) = 0.0;
end
if nyd==1,
	JC=ncy;
	KC=1:ncz;
	IC=1:ncx;
	rhsc(IC,JC,KC) = 0.0;
end
if nze==1,
	KC=1;
	JC=1:ncy;
	IC=1:ncx;
	rhsc(IC,JC,KC) = 0.0;
end
if nzf==1,
	KC=ncz;
	JC=1:ncy;
	IC=1:ncx;
	rhsc(IC,JC,KC) = 0.0;
end

function coarse=Full_Weighting_3d(fine,IM1,I,IP1,JM1,J,JP1,KM1,K,KP1, ...
	coarse,IC,JC,KC)

coarse(IC,JC,KC) = 0.25*( ...
	(fine(IM1,JM1,KM1)+fine(IP1,JM1,KM1)+fine(IM1,JP1,KM1)+fine(IP1,JP1,KM1)+...
	2.0*(fine(IM1,J,KM1)+fine(IP1,J,KM1)+fine(I,JM1,KM1)+fine(I,JP1,KM1))+...
	4.0*fine(I,J,KM1))*0.0625+ ...
	2.0*(fine(IM1,JM1,K)+fine(IP1,JM1,K)+fine(IM1,JP1,K)+fine(IP1,JP1,K)+...
	2.0*(fine(IM1,J,K)+fine(IP1,J,K)+fine(I,JM1,K)+fine(I,JP1,K))+...
	4.0*fine(I,J,K))*0.0625+ ...
	(fine(IM1,JM1,KP1)+fine(IP1,JM1,KP1)+fine(IM1,JP1,KP1)+fine(IP1,JP1,KP1)+...
	2.0*(fine(IM1,J,KP1)+fine(IP1,J,KP1)+fine(I,JM1,KP1)+fine(I,JP1,KP1))+...
	4.0*fine(I,J,KP1))*0.0625);

