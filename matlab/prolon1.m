function q=prolon1(p,q,k,solver)
% function q=prolon1(p,q,k,solver)

%
% $Id: prolon1.m,v 1.3 2011/03/26 12:56:39 patrick Exp $
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

[ncx,ncy,ncz]=size(p);
ncx=ncx-2;

xoff=solver.xoff;

nxa=solver.nxa;
nxb=solver.nxb;

intpol=solver.intpol;

nx=solver.nk{k}(1);

ist=1;
ioddst=1;
ifn=nx;
ioddfn=nx;
if nxa==1,
	ist=2; 
	ioddst=3; 
end;
if nxb==1, 
	ifn=nx-1; 
	ioddfn=nx-2; 
end;

if intpol==1 | ncx<4, % linearly interpolate in y

  if ncx<nx,
% every other point of nx grid is ncx grid
	  I=ioddst:2:ioddfn;
		IC=fix((I+1)/2);

		I=I+xoff; IC=IC+xoff;
		q(I,:,:)=p(IC,:,:);
		I=2:2:ifn;
		I=I+xoff;
		q(I,:,:)=.5*(q(I-1,:,:)+q(I+1,:,:));
	else
% nx grid equals ncx grid
		I=ist:ifn;
		I=I+xoff;
		q(I,:,:)=p(I,:,:);
	end

% set virtual end points if periodic
	if nxa==0,
		q(0+xoff,:,:)=q(nx-1+xoff,:,:);
		q(nx+1+xoff,:,:)=q(2+xoff,:,:);
	end
	return
else, % cubically interpolate in x

	if ncx<nx,

		I=ioddst:2:ioddfn;
		IC=fix((I+1)/2);
		I=I+xoff; IC=IC+xoff;
		q(I,:,:)=p(IC,:,:);

% set deep interior with symmetric formula

		I=4:2:nx-3;
		I=I+xoff;
		q(I,:,:)=(-q(I-3,:,:)+9.*(q(I-1,:,:)+q(I+1,:,:))-q(I+3,:,:))*.0625;

% interpolate from q at i=2 and i=nx-1
		if nxa ~= 0,
% asymmetric formula near nonperiodic bndys
			q(2+xoff,:,:)=(5.*q(1+xoff,:,:)+15.*q(3+xoff,:,:)-...
				5.0*q(5+xoff,:,:)+q(7+xoff,:,:))*.0625;
			q(nx-1+xoff,:,:)=(5.*q(nx+xoff,:,:)+15.*q(nx-2+xoff,:,:)-...
				5.*q(nx-4+xoff,:,:)+q(nx-6+xoff,:,:))*.0625;
		else
% periodicity in x alows symmetric formula near bndys
			q(2+xoff,:,:)=(-q(nx-2+xoff,:,:)+...
				9.*(q(1+xoff,:,:)+q(3+xoff,:,:))-q(5+xoff,:,:))*.0625;
			q(nx-1+xoff,:,:)=(-q(nx-4+xoff,:,:)+...
				9.*(q(nx-2+xoff,:,:)+q(nx+xoff,:,:))-q(3+xoff,:,:))*.0625;
			q(nx+1+xoff,:,:)=q(2+xoff,:,:);
			q(0+xoff,:,:)=q(nx-1+xoff,:,:);
		end
		return
	else
% ncx grid equals nx grid
		I=ist:ifn;
		q(I+xoff,:,:)=p(I+xoff,:,:);
		if nxa==0,
			q(0+xoff,:,:)=q(nx-1+xoff,:,:);
			q(nx+1+xoff,:,:)=q(2+xoff,:,:);
		end
		return
	end
end

