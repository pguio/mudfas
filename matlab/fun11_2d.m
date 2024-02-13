function [u,f,solver]=fun11(solver)
% function [u,f,solver]=fun11(solver)

%
% $Id: fun11_2d.m,v 1.2 2011/03/26 12:56:39 patrick Exp $
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

dust2file = '/home/patrick/research/publications/dust2/matlab/dust4_2dn6';
if exist([dust2 '.mat'],'file'),

eval(['load ' dust2file]);

solver.nx=size(rho,1);
solver.ny=size(rho,2);

solver.xa=-50;
solver.xb=50;
solver.yc=-50;
solver.yd=50;

solver.nxa=2;
solver.nxb=2;
solver.nyc=2;
solver.nyd=2;


solver_const;
solver.teshape=Te2dModel.Const2d;

solver.temin=1.0;
solver.temax=1.0;

solver.intpol=3;

u=zeros(solver.nx,solver.ny);

f=-double(rho);

else

error(['mat file ' dust2file ' not found.']);

end
