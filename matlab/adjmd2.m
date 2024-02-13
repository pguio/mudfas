function [phi,rhs]=adjmd2(phi,rhs,k,solver)
% function [phi,rhs]=adjmd2(phi,rhs,k,solver)
% adjust righ thand side for various boundary conditions

%
% $Id: adjmd2.m,v 1.5 2011/03/26 12:56:39 patrick Exp $
%
% Copyright (c) 2000-2011 Patrick Guio <patrick.guio@gmail.com>
% All Rights Reserved.
%
% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the Free
% Software Foundation; either version 2.  of the License, or (at your
% option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
% or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
% for more details.
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

xoff=solver.xoff; 
yoff=solver.yoff;

nx=solver.nk{k}(1);
ny=solver.nk{k}(2);

dlx=(xb-xa)/(nx-1);
dly=(yd-yc)/(ny-1);
dlx2=dlx+dlx;
dly2=dly+dly;
dlxx=dlx*dlx;
dlyy=dly*dly;

ist=1;
ifn=nx;
jst=1;
jfn=ny;
if nxa==1,
	ist=2; 
end
if nxb==1, 
	ifn=nx-1; 
end
if nyc==1, 
	jst=2; 
end
if nyd==1,
	jfn=ny-1; 
end

% adjust right hand side at derivative boundaries
if nxa==2,
  x=xa;
  I=1;
  cxx=1;
  cx=0;
  cxx=max(cxx,abs(cx)*dlx*0.5);
  c1=cxx/dlxx-cx/dlx2;
  J=jst:jfn;
  Y=yc+(J-1)*dly;
  gbdy=solver.gbdxa;
	rhs(I,J)=rhs(I,J)+dlx2*c1*gbdy;
end
if nxb==2,
  x=xb;
  I=nx;
  cxx=1;
  cx=0;
  cxx=max(cxx,abs(cx)*dlx*0.5);
  c2=cxx/dlxx+cx/dlx2;
  J=jst:jfn;
  Y=yc+(J-1)*dly;
  gbdy=solver.gbdxb;
  rhs(I,J)=rhs(I,J)-dlx2*c2*gbdy;
end
if nyc==2,
  y=yc;
  J=1;
  cyy=1;
  cy=0;
  cyy=max(cyy,abs(cy)*dly*0.5);
  c1=cyy/dlyy-cy/dly2;
  I=ist:ifn;
  X=xa+(I-1)*dlx;
  gbdy=solver.gbdyc;
  rhs(I,J)=rhs(I,J)+dly2*c1*gbdy;
end
if nyd==2,
  y=yd;
  J=ny;
  cyy=1;
  cy=0;
  cyy=max(cyy,abs(cy)*dly*0.5);
  c2=cyy/dlyy+cy/dly2;
  I=ist:ifn;
  X=xa+(I-1)*dlx;
  gbdy=solver.gbdyd;
  rhs(I,J)=rhs(I,J)-dly2*c2*gbdy;
end

% set specified boundaries in rhs from phi
if nxa==1,
	I=1;
	J=1:ny;
	rhs(I,J)=phi(I+xoff,J+yoff);
end
if nxb==1,
	I=nx;
	J=1:ny;
	rhs(I,J)=phi(I+xoff,J+yoff);
end
if nyc==1,
	J=1;
	I=1:nx;
	rhs(I,J)=phi(I+xoff,J+yoff);
end
if nyd==1,
	J=ny;
	I=1:nx;
	rhs(I,J)=phi(I+xoff,J+yoff);
end
