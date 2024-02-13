function [phi,rhs]=adjmd3(phi,rhs,k,solver)
% function [phi,rhs]=adjmd3(phi,rhs,k,solver)
% adjust righthand side for various boundary conditions

%
% $Id: adjmd3.m,v 1.4 2011/03/26 12:56:39 patrick Exp $
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
ze=solver.ze;
zf=solver.zf;

nxa=solver.nxa;
nxb=solver.nxb;
nyc=solver.nyc;
nyd=solver.nyd;
nze=solver.nze;
nzf=solver.nzf;

xoff=solver.xoff;
yoff=solver.yoff;
zoff=solver.zoff;

nx=solver.nk{k}(1);
ny=solver.nk{k}(2);
nz=solver.nk{k}(3);

dlx=(xb-xa)/(nx-1);
dlx2=dlx+dlx;
dlxx=dlx*dlx;
dly=(yd-yc)/(ny-1);
dly2=dly+dly;
dlyy=dly*dly;
dlz=(zf-ze)/(nz-1);
dlz2=dlz+dlz;
dlzz=dlz*dlz;

ist=1;
ifn=nx;
jst=1;
jfn=ny;
kst=1;
kfn=nz;
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
if nze==1, 
	kst=2; 
end
if nzf==1,
	kfn=nz-1; 
end

% adjust right hand side at derivative boundaries
if nxa==2,
	x=xa;
	I=1;
	cxx=1;
	cx=0;
	cxx=max(cxx,abs(cx)*dlx*0.5);
	c1=cxx/dlxx-cx/dlx2;
	K=kst:kfn;
	Z=ze+(K-1)*dlz;
	J=jst:jfn;
	Y=yc+(J-1)*dly;
	gbdy=solver.gbdxa;
	rhs(I,J,K)=rhs(I,J,K)+dlx2*c1*gbdy;
end
if nxb==2,
	x=xb;
	I=nx;
	cxx=1;
	cx=0;
	cxx=max(cxx,abs(cx)*dlx*0.5);
	c2=cxx/dlxx+cx/dlx2;
	K=kst:kfn;
	Z=ze+(K-1)*dlz;
	J=jst:jfn;
	Y=yc+(J-1)*dly;
	gbdy=solver.gbdxb;
	rhs(I,J,K)=rhs(I,J,K)-dlx2*c2*gbdy;
end
if nyc==2,
	y=yc;
	J=1;
	cyy=1;
	cy=0;
	cyy=max(cyy,abs(cy)*dly*0.5);
	c1=cyy/dlyy-cy/dly2;
	K=kst:kfn;
	Z=ze+(K-1)*dlz;
	I=ist:ifn;
	X=xa+(I-1)*dlx;
	gbdy=solver.gbdyc;
	rhs(I,J,K)=rhs(I,J,K)+dly2*c1*gbdy;
end
if nyd==2,
	y=yd;
	J=ny;
	cyy=1;
	cy=0;
	cyy=max(cyy,abs(cy)*dly*0.5);
	c2=cyy/dlyy+cy/dly2;
	K=kst:kfn;
	Z=ze+(K-1)*dlz;
	I=ist:ifn;
	X=xa+(I-1)*dlx;
	gbdy=solver.gbdyd;
	rhs(I,J,K)=rhs(I,J,K)-dly2*c2*gbdy;
end
if nze==2,
	z=ze;
	K=1;
	czz=1;
	cz=0;
	czz=max(czz,abs(cz)*dlz*0.5);
	c1=czz/dlzz-cz/dlz2;
	J=jst:jfn;
	Y=yc+(J-1)*dly;
	I=ist:ifn;
	X=xa+(I-1)*dlx;
	gbdy=solver.gbdze;
	rhs(I,J,K)=rhs(I,J,K)+dlz2*c1*gbdy;
end
if nzf==2,
	z=zf;
	K=nz;
	czz=1;
	cz=0;
	czz=max(czz,abs(cz)*dlz*0.5);
	c2=czz/dlzz+cz/dlz2;
	J=jst:jfn;
	Y=yc+(J-1)*dlz;
	I=ist:ifn;
	X=xa+(I-1)*dlx;
	gbdy=solver.gbdzf;
	rhs(I,J,K)=rhs(I,J,K)-dlz2*c2*gbdy;
end

% set specified boundaries in rhs from phi
if nxa==1,
	I=1;
	J=1:ny;
	K=1:nz;
	rhs(I,J,K)=phi(I+xoff,J+yoff,K+zoff);
end
if nxb==1,
	I=nx;
	J=1:ny;
	K=1:nz;
	rhs(I,J,K)=phi(I+xoff,J+yoff,K+zoff);
end
if nyc==1,
	J=1;
	K=1:nz;
	I=1:nx;
	rhs(I,J,K)=phi(I+xoff,J+yoff,K+zoff);
end
if nyd==1,
	J=ny;
	K=1:nz;
	I=1:nx;
	rhs(I,J,K)=phi(I+xoff,J+yoff,K+zoff);
end
if nze==1,
	K=1;
	J=1:ny;
	I=1:nx;
	rhs(I,J,K)=phi(I+xoff,J+yoff,K+zoff);
end
if nzf==1,
	K=nz;
	J=1:ny;
	I=1:nx;
	rhs(I,J,K)=phi(I+xoff,J+yoff,K+zoff);
end

