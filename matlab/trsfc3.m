function [phic,rhsc]=trsfc3(phi,rhs,kc,solver)
% function [phic,rhsc]=trsfc3(phi,rhs,kc,solver)
% transfer fine grid to coarse grid

%
% $Id: trsfc3.m,v 1.3 2011/03/26 12:56:40 patrick Exp $
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

nx=solver.nk{kc+1}(1);
ny=solver.nk{kc+1}(2);
nz=solver.nk{kc+1}(3);

ncx=solver.nk{kc}(1);
ncy=solver.nk{kc}(2);
ncz=solver.nk{kc}(3);

xoff=solver.xoff;
yoff=solver.yoff;
zoff=solver.zoff;

% initialise coarse phic
% and set virtual boundaries in phic to zero
phic=zeros(ncx+2,ncy+2,ncz+2);
% initialise coarse rhs
rhsc=zeros(ncx,ncy,ncz);

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

IC=1:ncx;
I=IC+ix*(IC-1);

JC=1:ncy;
J=JC+jy*(JC-1);

KC=1:ncz;
K=KC+kz*(KC-1);

phic(IC+xoff,JC+yoff,KC+zoff)=phi(I+xoff,J+yoff,K+zoff);
rhsc(IC,JC,KC)=rhs(I,J,K);

