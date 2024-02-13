function solver=initpde2(solver)
% function solver=initpde2(solver)

%
% $Id: initpde2.m,v 1.3 2011/03/26 12:56:39 patrick Exp $
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

[ixp,iex]=calcnmudpack(solver.nx);
[jyq,jey]=calcnmudpack(solver.ny);

solver.ixp=ixp;
solver.jyq=jyq; 
solver.iex=iex; 
solver.jey=jey; 
solver.ngrid=max([iex jey]); 

ngrid=solver.ngrid;

iex=solver.iex;
ixp=solver.ixp;
jey=solver.jey;
jyq=solver.jyq;

xa=solver.xa;
xb=solver.xb;
yc=solver.yc;
yd=solver.yd;

for k=1:ngrid,
% Calculate grid size at each level
  nx = ixp*2^(max(k+iex-ngrid,1)-1)+1;
	ny = jyq*2^(max(k+jey-ngrid,1)-1)+1;

	solver.nk{k}=[nx,ny];

% Damping parameterof coarse grid-correction
	solver.sk{k}=1.0;
end

nxa=solver.nxa;
nxb=solver.nxb;
nyc=solver.nyc;
nyd=solver.nyd;

for k=1:ngrid,

	nx=solver.nk{k}(1);
	ny=solver.nk{k}(2);
	dlx=(xb-xa)/(nx-1);
	dlx2=dlx+dlx;
	dlxx=dlx*dlx;
	dly=(yd-yc)/(ny-1);
	dly2=dly+dly;
	dlyy=dly*dly;
	cmin=1.0;
	alfmax=0.0;
	cemax=0.0;

% set x,y subscript limits for calls to cofx,cofy
% (avoid specified boundaries)
	ist=1;
	ifn=nx;
	jst=1;
	jfn=ny;
	if nxa==1,
		ist=2;
	end;
	if nxb==1,
		ifn=nx-1;
	end;
	if nyc==1,
		jst=2;
	end;
	if nyd==1,
		jfn=ny-1;
	end;

% compute discretization coefficients on interior and
% nonspecified boundaries
	J=jst:jfn;
	Y=yc+(J-1)*dly;
	cyy=1.0;
	cy=0.0;
	cey=0.0;
	cmin=min(cmin,cyy);
	cemax=max(abs(cey),cemax);

	cyy=max(cyy,abs(cy)*dly*0.5);
	c1=cyy/dlyy-cy/dly2;
	c2=cyy/dlyy+cy/dly2;
	c3=cey-(c1+c2);
	cofy(J,1)=c1;
	cofy(J,2)=c2;
	cofy(J,3)=c3;

	I=ist:ifn;
	X=xa+(I-1)*dlx;
	cxx=1.0;
	cx=0.0;
	cex=0.0;
	cmin=min(cmin,cxx);
	cemax=max(abs(cex),cemax);

	cxx=max(cxx,abs(cx)*dlx*0.5);
	c1=cxx/dlxx-cx/dlx2;
	c2=cxx/dlxx+cx/dlx2;
	c3=cex-(c1+c2);
	cofx(I,1)=c1;
	cofx(I,2)=c2;
	cofx(I,3)=c3;

% adjust discretization for mixed derivative b.c.
	if nxa==2,
		I=1;
		c1=cofx(I,1);
		cofx(I,1)=0.0;
		cofx(I,2)=cofx(I,2)+c1;
		Y=yc+dly;
% compute constant coefficient alfa
		alfa=solver.alfxa;
		alfmax=max(alfmax,abs(alfa));
		cofx(I,3)=cofx(I,3)+dlx2*alfa*c1;
	end
	if nxb==2,
		I=nx;
		Y=yc+dly;
% compute constant coefficient alfa
		alfa=solver.alfxb;
		c2=cofx(I,2);
		cofx(I,1)=cofx(I,1)+c2;
		cofx(I,2)=0.0;
		cofx(I,3)=cofx(I,3)-dlx2*alfa*c2;
		alfmax=max(alfmax,abs(alfa));
	end
	if nyc==2,
		J=1;
		X=xa+dlx;
% compute constant coefficient alfa
		alfa=solver.alfyc;
		c1=cofy(J,1);
		cofy(J,1)=0.0;
		cofy(J,2)=cofy(J,2)+c1;
		cofy(J,3)=cofy(J,3)+dly2*alfa*c1;
		alfmax=max(alfmax,abs(alfa));
	end
	if nyd==2,
		J=ny;
		X=xa+dlx;
% compute constant coefficient alfa
		alfa=solver.alfyd;
		c2=cofy(J,2);
		cofy(J,2)=0.0;
		cofy(J,1)=cofy(J,1)+c2;
		cofy(J,3)=cofy(J,3)-dly2*alfa*c2;
		alfmax=max(alfmax,abs(alfa));
	end

	solver.cofx{k}=cofx;
	solver.cofy{k}=cofy;

end

