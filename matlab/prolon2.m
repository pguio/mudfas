function q=prolon2(p,q,k,solver)
% function q=prolon2(p,q,k,solver)

%
% $Id: prolon2.m,v 1.3 2011/03/26 12:56:39 patrick Exp $
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

[ncx,ncy,nz]=size(p);
ncx=ncx-2;
ncy=ncy-2;

xoff=solver.xoff; 
yoff=solver.yoff;

nxa=solver.nxa;
nxb=solver.nxb;
nyc=solver.nyc;
nyd=solver.nyd;

intpol=solver.intpol;

nx=solver.nk{k}(1);
ny=solver.nk{k}(2);

ist=1;
ifn=nx;
jst=1;
jfn=ny;
joddst=1;
joddfn=ny;
if nxa==1, 
	ist=2; 
end;
if nxb==1, 
	ifn=nx-1; 
end;
if nyc==1,
	jst=2; 
	joddst=3; 
end;
if nyd==1,
	jfn=ny-1;
	joddfn=ny-2; 
end;

if intpol==1 | ncy<4, % linearly interpolate in y

	if ncy<ny,
% ncy grid is an every other point subset of ny grid
% set odd j lines interpolating in x and then set even
% j lines by averaging odd j lines
		J=joddst:2:joddfn;
		JC=fix(J/2+1);
		J=J+yoff; JC=JC+yoff;
		q(:,J,:)=prolon1(p(:,JC,:),q(:,J,:),k,solver);
		J=2:2:jfn;
		I=ist:ifn;
		J=J+yoff; I=I+xoff;
		q(I,J,:)=.5*(q(I,J-1,:)+q(I,J+1,:));

% set periodic virtual boundaries if necessary
    if nyc==0,
		  I=ist:ifn;
			q(I+xoff,0+yoff,:)=q(I+xoff,ny-1+yoff,:);
			q(I+xoff,ny+1+yoff,:)=q(I+xoff,2+yoff,:);
		end
		return
	else,
% ncy grid equals ny grid so interpolate in x only
		J=jst:jfn;
		JC=J;
		J=J+yoff; JC=JC+yoff;
		q(:,J,:)=prolon1(p(:,JC,:),q(:,J,:),k,solver);

% set periodic virtual boundaries if necessary
    if nyc==0,
		  I=ist:ifn;
			I=I+xoff;
			q(I,0+yoff,:)=q(I,ny-1+yoff,:);
			q(I,ny+1+yoff,:)=q(I,2+yoff,:);
		end
		return
	end
else, % cubically interpolate in y

	if ncy<ny,
% set every other point of ny grid by interpolating in x
		J=joddst:2:joddfn;
		JC=fix(J/2+1);
		J=J+yoff; JC=JC+yoff;
		q(:,J,:)=prolon1(p(:,JC,:),q(:,J,:),k,solver);

% set deep interior of ny grid using values just
% generated and symmetric cubic interpolation in y

		J=4:2:ny-3;
		I=ist:ifn;
		I=I+xoff; J=J+yoff;
		q(I,J,:)=(-q(I,J-3,:)+9.*(q(I,J-1,:)+q(I,J+1,:))-q(I,J+3,:))*.0625;

% interpolate from q at j=2 and j=ny-1
    if nyc ~= 0,
% asymmetric formula near nonperiodic y boundaries
			I=ist:ifn;
			I=I+xoff;
			q(I,2+yoff,:)=(5.*q(I,1+yoff,:)+15.0*q(I,3+yoff,:)- ...
				5.*q(I,5+yoff,:)+q(I,7+yoff,:))*.0625;
			q(I,ny-1+yoff,:)=(5.*q(I,ny+yoff,:)+15.*q(I,ny-2+yoff,:)- ...
				5.*q(I,ny-4+yoff,:)+q(I,ny-6+yoff,:))*.0625;
		else
% periodicity in y alows symmetric formula near bndys
			I=ist:ifn;
			I=I+xoff;
			q(I,2+yoff,:)=(-q(I,ny-2+yoff,:)+9.*(q(I,1+yoff,:)+ ...
				q(I,3+yoff,:))-q(I,5+yoff,:))*.0625;
			q(I,ny-1+yoff,:)=(-q(I,ny-4+yoff,:)+9.*(q(I,ny-2+yoff,:)+ ...
				q(I,ny+yoff,:))-q(I,3+yoff,:))*.0625;
			q(I,ny+1+yoff,:)=q(I,2+yoff,:);
			q(I,0+yoff,:)=q(I,ny-1+yoff,:);
		end
		return
	else
% ncy grid is equals ny grid so interpolate in x only
		J=jst:jfn;
		JC=J;
		J=J+yoff; JC=JC+yoff;
		q(:,J,:)=prolon1(p(:,JC,:),q(:,J,:),k,solver);
% set periodic virtual boundaries if necessary
		if nyc==0,
			I=ist:ifn;
			I=I+xoff;
			q(I,0+yoff,:)=q(I,ny-1+yoff,:);
			q(I,ny+1+yoff,:)=q(I,2+yoff,:);
		end
		return
	end
end

