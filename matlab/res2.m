function rhsc=res2(resf,rhsc,k,solver)
% function rhsc=res2(resf,rhsc,k,solver)
%
% restrict resf grid residual in resf to coarse grid in rhsc
% using full weighting for all double precision 2d codes

%
% $Id: res2.m,v 1.3 2011/03/26 12:56:39 patrick Exp $
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

nx=solver.nk{k+1}(1);
ny=solver.nk{k+1}(2);

ncx=solver.nk{k}(1);
ncy=solver.nk{k}(2);

% set x,y coarsening integer subscript scales
ix=1;
if ncx==nx,
	ix=0; 
end;
jy=1;
if ncy==ny,
	jy=0;
end;

% restrict on interior
IC=2:ncx-1;
I=IC+ix*(IC-1);
JC=2:ncy-1;
J=JC+jy*(JC-1);

IM1=I-1; IP1=I+1;
JM1=J-1; JP1=J+1;
rhsc=Full_Weighting_2d(resf,IM1,I,IP1,JM1,J,JP1,rhsc,IC,JC);

% set residual on boundaries

% y=yc,yd boundaries
JC=[1 ncy];
J=JC+jy*(JC-1);
JM1=max(J-1,2);
JP1=min(J+1,ny-1);
JM1(find(J==1 & nyc==0))=ny-1;
JP1(find(J==ny & nyc==0))=2;

% y=yc,yd and x=xa,xb corners
IC=[1 ncx];
I=IC+ix*(IC-1);
IM1=max(I-1,2);
IP1=min(I+1,nx-1);
IM1(find(I==1 & nxa==0))=nx-1;
IP1(find(I==nx & nxb==0))=2;

rhsc=Full_Weighting_2d(resf,IM1,I,IP1,JM1,J,JP1,rhsc,IC,JC);

% set y=yc,yd interior edges
IC=2:ncx-1;
I=IC+ix*(IC-1);

IM1=I-1; IP1=I+1;
rhsc=Full_Weighting_2d(resf,IM1,I,IP1,JM1,J,JP1,rhsc,IC,JC);

% set x=xa,xb interior edges
IC=[1 ncx];
I=IC+ix*(IC-1);
IM1 = max(I-1,2);
IP1 = min(I+1,nx-1);
IM1(find(I==1 & nxa==0))=nx-1;
IP1(find(I==nx & nxb==0))=2;

JC=2:ncy-1;
J=JC+jy*(JC-1);

JM1=J-1; JP1=J+1;
rhsc=Full_Weighting_2d(resf,IM1,I,IP1,JM1,J,JP1,rhsc,IC,JC);

return;

% set coarse grid residual zero on specified boundaries
if nxa==1,
	IC=1;
	JC=1:ncy;
  rhsc(IC,JC)=0.0;
end
if nxb==1,
	IC=ncx;
	JC=1:ncy;
  rhsc(IC,JC)=0.0;
end
if nyc==1,
	IC=1:ncx;
	JC=1;
  rhsc(IC,JC)=0.0;
end
if nyd==1,
	IC=1:ncx;
	JC=ncy;
  rhsc(IC,JC)=0.0;
end

function coarse=Full_Weighting_2d(fine,IM1,I,IP1,JM1,J,JP1, ...
	coarse, IC, JC)

coarse(IC,JC)=( ...
	fine(IM1,JM1)+fine(IP1,JM1)+fine(IM1,JP1)+fine(IP1,JP1)+ ...
	2.0*(fine(IM1,J)+fine(IP1,J)+fine(I,JM1)+fine(I,JP1))+ ...
	4.0*fine(I,J))*0.0625;


