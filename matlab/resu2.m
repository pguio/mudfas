function phic=resu2(phif,phic,k,solver)
% function phic=resu2(phif,phic,k,solver)
%
% restrict phi grid residual in phif to coarse grid in phic
% using full weighting for all double precision 2d codes

%
% $Id: resu2.m,v 1.3 2011/03/26 12:56:39 patrick Exp $
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

xoff=solver.xoff;
yoff=solver.yoff;

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

JCo=JC+yoff; Jo=J+yoff;
ICo=IC+xoff; Io=I+xoff;

IM1o=Io-1; IP1o=Io+1;
JM1o=Jo-1; JP1o=Jo+1;
phic=Full_Weighting_2d(phif,IM1o,Io,IP1o,JM1o,Jo,JP1o,phic,ICo,JCo);

% set residual on boundaries

% y=yc,yd boundaries
JC=[1 ncy];
J=JC+jy*(JC-1);
JM1=max(J-1,2);
JP1=min(J+1,ny-1);
JM1(find(J==1 & nyc==0))=ny-1;
JP1(find(J==ny & nyd==0))=2;

% y=yc,yd and x=xa,xb corners
IC=[1 ncx];
I=IC+ix*(IC-1);
IM1=max(I-1,2);
IP1=min(I+1,nx-1);
IM1(find(I==1 & nxa==0))=nx-1;
IP1(find(I==nx & nxb==0))=2;

JCo=JC+yoff; Jo=J+yoff; JM1o=JM1+yoff; JP1o=JP1+yoff;
ICo=IC+xoff; Io=I+xoff; IM1o=IM1+xoff; IP1o=IP1+xoff;

phic=Full_Weighting_2d(phif,IM1o,Io,IP1o,JM1o,Jo,JP1o,phic,ICo,JCo);

% set y=yc,yd interior edges
IC=2:ncx-1;
I=IC+ix*(IC-1);

JCo=JC+yoff; Jo=J+yoff; JM1o=JM1+yoff; JP1o=JP1+yoff;
ICo=IC+xoff; Io=I+xoff;

IM1o=Io-1; IP1o=Io+1;
phic=Full_Weighting_2d(phif,IM1o,Io,IP1o,JM1o,Jo,JP1o,phic,ICo,JCo);

% set x=xa,xb interior edges
IC=[1 ncx];
I=IC+ix*(IC-1);
IM1 = max(I-1,2);
IP1 = min(I+1,nx-1);
IM1(find(I==1 & nxa==0))=nx-1;
IP1(find(I==nx & nxb==0))=2;
JC=2:ncy-1;
J=JC+jy*(JC-1);

JCo=JC+yoff; Jo=J+yoff;
ICo=IC+xoff; Io=I+xoff; IM1o=IM1+xoff; IP1o=IP1+xoff;

JM1o=Jo-1; JP1o=Jo+1;
phic=Full_Weighting_2d(phif,IM1o,Io,IP1o,JM1o,Jo,JP1o,phic,ICo,JCo);

return;

% set coarse grid residual zero on specified boundaries
if nxa==1,
  IC=1;
	JC=1:ncy;
  phic(IC+xoff,JC+yoff)=0.0;
end
if nxb==1,
	IC=ncx;
	JC=1:ncy;
  phic(IC+xoff,JC+yoff)=0.0;
end
if nyc==1,
	IC=1:ncx;
	JC=1;
  phic(IC+xoff,JC+yoff)=0.0;
end
if nyd==1,
	IC=1:ncx;
	JC=ncy;
  phic(IC+xoff,JC+yoff)=0.0;
end

function coarse=Full_Weighting_2d(fine,IM1,I,IP1,JM1,J,JP1, ...
	coarse, IC, JC)

coarse(IC,JC)=( ...
	fine(IM1,JM1)+fine(IP1,JM1)+fine(IM1,JP1)+fine(IP1,JP1)+ ...
	2.0*(fine(IM1,J)+fine(IP1,J)+fine(I,JM1)+fine(I,JP1))+ ...
	4.0*fine(I,J))*0.0625;


