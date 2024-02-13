function [phi,rhs]=swk3(phif,rhsf,solver)
% function [phi,rhs]=swk3(phif,rhsf,solver)
%
% set phif,rhsf input in arrays which include
% virtual boundaries for phi (for all 2-d double precision codes)

%
% $Id: swk3.m,v 1.3 2011/03/26 12:56:39 patrick Exp $
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

nfx=solver.nk{solver.ngrid}(1);
nfy=solver.nk{solver.ngrid}(2);
nfz=solver.nk{solver.ngrid}(3);

xoff=solver.xoff; 
yoff=solver.yoff; 
zoff=solver.zoff;

% initialise phi
% and set virtual boundaries in phi to zero
phi=zeros(nfx+2,nfy+2,nfz+2);

I=1:nfx;
J=1:nfy;
K=1:nfz;

phi(I+xoff,J+yoff,K+zoff)=phif(I,J,K);

% change sign of ni in order to solve
% \nabla^2 phi - exp(u/Te) = -ni
rhs=zeros(nfx,nfy,nfz);
rhs(I,J,K)=-rhsf(I,J,K);


