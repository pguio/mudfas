function [u,f,solver]=mgfas2d(f,varargin)
% function [u,f,solver]=mgfas2d(f,varargin)

%
% $Id: mgfas2d.m,v 1.4 2011/03/26 12:56:39 patrick Exp $
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

solver=solver2d_def;
solver=solver_par_parsing(solver,varargin{:});

[u,f,solver]=feval([f '_2d'],solver);

% discretize pde at each grid level
solver=initpde2(solver);

view_solver(solver);

lim=[solver.xa solver.xb solver.yc solver.yd];
te=[solver.temin,solver.temax,solver.sigmat];

mgfasflag=0;
argc=length(varargin);
for i=1:2:argc,
	param=varargin{i};
	val=varargin{i+1};
	switch lower(param)
		case 'mgfasflag', mgfasflag=val;
	end
end

[u,te]=MGFAS2D(f,lim,te,mgfasflag);

f=-f;

plotdefect2d(u,f,solver)
