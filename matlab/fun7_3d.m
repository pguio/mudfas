function [u,f,solver]=fun9(solver)
% function [u,f,solver]=fun9(solver)

%
% $Id: fun7_3d.m,v 1.2 2011/03/26 12:56:39 patrick Exp $
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

solver.nx=57;
solver.ny=57;
solver.nz=57;

if 1,
solver.nxa=1;
solver.nxb=1;
solver.nyc=1;
solver.nyd=1;
solver.nze=1;
solver.nzf=1;
end


% Initialise phi
u=zeros(solver.nx,solver.ny,solver.nz);

% Initialise rho
f=ones(solver.nx,solver.ny,solver.nz);

f(fix(solver.nx/2), fix(solver.ny/2), fix(solver.nz/2)) = 300;
