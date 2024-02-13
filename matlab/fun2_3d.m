function [u,f,solver]=fun2(solver)
% function [u,f,solver]=fun2(solver)

%
% $Id: fun2_3d.m,v 1.2 2011/03/26 12:56:39 patrick Exp $
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

u=zeros(solver.nx,solver.ny,solver.nz);

f=ones(solver.nx,solver.ny,solver.nz);

imn=fix(solver.nx/3); 
imx=fix(2*solver.nx/3);
jmn=fix(solver.ny/3); 
jmx=fix(2*solver.ny/3);
kmn=fix(solver.nz/3); 
kmx=solver.nz;

f(imn:imx,jmn:jmx,kmn:kmx)=1.2*f(imn:imx,jmn:jmx,kmn:kmx);


