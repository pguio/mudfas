function q=prolon3(p,q,k,solver)
% function q=prolon3(p,q,k,solver)

%
% $Id: prolon3.m,v 1.3 2011/03/26 12:56:39 patrick Exp $
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

ncx=solver.nk{k-1}(1);
ncy=solver.nk{k-1}(2);
ncz=solver.nk{k-1}(3);

xoff=solver.xoff;
yoff=solver.yoff;
zoff=solver.zoff;

nxa=solver.nxa;
nxb=solver.nxb;
nyc=solver.nyc;
nyd=solver.nyd;
nze=solver.nze;
nzf=solver.nzf;

intpol=solver.intpol;

nx=solver.nk{k}(1);
ny=solver.nk{k}(2);
nz=solver.nk{k}(3);

ist=1;
ifn=nx;
jst=1;
jfn=ny;
kst=1;
kfn=nz;
koddst=1;
koddfn=nz;
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
if nze==1,
	kst=2;
	koddst=3;
end;
if nzf==1,
	kfn=nz-1;
	koddfn=nz-2;
end

if intpol==1 | ncz<4, 
% linearly interpolate in z

  if ncz<nz,
% ncz grid is an every other point subset of nz grid
% set odd k planes interpolating in x&y and then set even
% k planes by averaging odd k planes
		K=koddst:2:koddfn;
		KC=fix(K/2+1);
		K=K+zoff; KC=KC+zoff;
		q(:,:,K)=prolon2(p(:,:,KC),q(:,:,K),k,solver);
		K=2:2:kfn;
		J=jst:jfn;
		I=ist:ifn;
		K=K+zoff; J=J+yoff; I=I+xoff;
		q(I,J,K)=.5*(q(I,J,K-1)+q(I,J,K+1));

% set periodic virtual boundaries if necessary
    if nze==0,
		  J=jst:jfn;
		  I=ist:ifn;
			J=J+yoff; I=I+xoff;
			q(I,J,0+zoff)=q(I,J,nz-1+zoff);
			q(I,J,nz+1+zoff)=q(I,J,2+zoff);
		end
		return
	else,
% ncz grid is equals nz grid so interpolate in x&y only
		K=kst:kfn;
		KC=K;
		K=K+zoff; KC=KC+zoff;
		q(:,:,K)=prolon2(p(:,:,KC),q(:,:,K),k,solver);

% set periodic virtual boundaries if necessary
    if nze==0,
		  J=jst:jfn;
		  I=ist:ifn;
			J=J+yoff; I=I+xoff;
			q(I,J,0+zoff)=q(I,J,nz-1+zoff);
			q(I,J,nz+1+zoff)=q(I,J,2+zoff);
		end
		return
	end
else, % cubically interpolate in z

	if ncz<nz,
% set every other point of nz grid by interpolating in x&y
		K=koddst:2:koddfn;
		KC=fix(K/2+1);
		K=K+zoff; KC=KC+zoff;
		q(:,:,K)=prolon2(p(:,:,KC),q(:,:,K),k,solver);

% set deep interior of nz grid using values just
% generated and symmetric cubic interpolation in z

		K=4:2:nz-3;
		J=jst:jfn;
		I=ist:ifn;
		J=J+yoff; I=I+xoff;
		q(I,J,K+zoff)=(-q(I,J,K+zoff-3)+9.*(q(I,J,K-1+zoff)+ ...
			q(I,J,K+1+zoff))-q(I,J,K+3+zoff))*.0625;

% interpolate from q at k=2 and k=nz-1
    if nze~=0,
% asymmetric formula near nonperiodic z boundaries
			J=jst:jfn;
			I=ist:ifn;
			J=J+yoff; I=I+xoff;
			q(I,J,2+zoff)=(5.*q(I,J,1+zoff)+15.*q(I,J,3+zoff)-5.*q(I,J,5+zoff)+ ...
				q(I,J,7+zoff))*.0625;
			q(I,J,nz-1+zoff)=(5.*q(I,J,nz+zoff)+15.*q(I,J,nz-2+zoff)- ...
				5.*q(I,J,nz-4+zoff)+q(I,J,nz-6+zoff))*.0625;
		else
% periodicity in y alows symmetric formula near bndys
			J=jst:jfn;
			I=ist:ifn;
			J=J+yoff; I=I+xoff;
			q(I,J,2+zoff)=(-q(I,J,nz-2+zoff)+9.*(q(I,J,1+zoff)+q(I,J,3+zoff))- ...
				q(I,J,5+zoff))*.0625;
			q(I,J,nz-1+zoff)=(-q(I,J,nz-4+zoff)+9.*(q(I,J,nz-2+zoff)+ ...
				q(I,J,nz+zoff))-q(I,J,3+zoff))*.0625;
			q(I,J,nz+1+zoff)=q(I,J,2+zoff);
			q(I,J,0+zoff)=q(I,J,nz-1+zoff);
		end
		return
	else

% ncz grid is equals nx grid so interpolate in x&y only
		K=kst:kfn;
		KC=K;
		K=K+xoff; KC=KC+xoff;
		q(:,:,K)=prolon2(p(:,:,K),q(:,:,K),k,solver);
% set periodic virtual boundaries if necessary
		if nze==0,
			J=jst:jfn;
			I=ist:ifn;
			J=J+yoff; I=I+xoff;
			q(I,J,0+zoff)=q(I,J,nz-1+zoff);
			q(I,J,nz+1+zoff)=q(I,J,2+zoff);
		end
		return
	end
end

