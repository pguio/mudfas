function [phic,rhsc]=trsfc2(phi,rhs,kc,solver)
% function [phic,rhsc]=trsfc2(phi,rhs,kc,solver)
%
% transfer fine grid to coarse grid

%
% $Id: trsfc2.m,v 1.3 2011/03/26 12:56:39 patrick Exp $
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

ncx=solver.nk{kc}(1);
ncy=solver.nk{kc}(2);

xoff=solver.xoff; 
yoff=solver.yoff;

% initialise coarse phic
% and set virtual boundaries in phic to zero
phic=zeros(ncx+2,ncy+2);
% initialise coarse rhs
rhsc=zeros(ncx,ncy);

ix=1;
if ncx==nx,
	ix=0;
end;

jy=1;
if ncy==ny,
	jy=0;
end;

IC=1:ncx;
I=IC+ix*(IC-1);

JC=1:ncy;
J=JC+jy*(JC-1);

phic(IC+xoff,JC+yoff)=phi(I+xoff,J+yoff);
rhsc(IC,JC)=rhs(I,J);

