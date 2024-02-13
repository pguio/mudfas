function solver=solver3d_def()
% function solver=solver3d_def()

%
% $Id: solver3d_def.m,v 1.5 2011/03/26 12:56:39 patrick Exp $
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


solver_const;

solver=struct('dim',[],'nx',[],'ny',[],'nz',[], ...
'nxa',[],'nxb',[],'nyc',[],'nyd',[],'nze',[],'nzf',[], ...
'xa',[],'xb',[],'yc',[],'yd',[],'ze',[],'zf',[], ...
'ixp',[],'jyq',[],'kzr',[],'iex',[],'jey',[],'kez',[],'ngrid',[], ...
'tolmax',[],'delta',[],'maxcy',[],'iguess',[], ...
'kcycle',[],'iprer',[],'ipost',[],'intpol',[], ...
'xoff',[],'yoff',[],'zoff',[], ...
'temin',[],'temax',[],'sigmat',[],...
'display',[],...
'nk',[],'bk',[],'sk',[]);

solver.dim=3;

solver.nx=37; 
solver.ny=37;
solver.nz=37;

% nxa,nxb,nyc,nyd is 
% 0 periodic
% 1 specified at boundary
% 2 mixed derivative boundary condition
solver.nxa=bc.mixed;
solver.nxb=bc.mixed;
solver.nyc=bc.mixed;
solver.nyd=bc.mixed;
solver.nze=bc.mixed;
solver.nzf=bc.mixed;

% domain limits
solver.xa=0;
solver.xb=36;
solver.yc=0;
solver.yd=36;
solver.ze=0;
solver.zf=36;

% Parameters for mixed boundary condition 
% For example at x=xa
% d(phi(xa, y, z))/dx + alfxa* phi(xa, y, z) = gbdxa(y, z)
solver.gbdxa=0;
solver.gbdxb=0;
solver.gbdyc=0;
solver.gbdyd=0;
solver.gbdze=0;
solver.gbdzf=0;
%
solver.alfxa=0;
solver.alfxb=0;
solver.alfyc=0;
solver.alfyd=0;
solver.alfze=0;
solver.alfzf=0;

% convergence paramters
solver.tolmax=5e-4;
solver.delta=5e-1;
solver.maxcy=10;
solver.iguess=0;

% kcycle = 1 corresponds to v cycles
% kcycle = 2 corresponds to w cycles
solver.kcycle=schedule.W;
solver.iprer=2;
solver.ipost=1;
% intpol = 1 corresponds to linear interpolation
% intpol = 3 corresponds to cubic interpolation
solver.intpol=interp.Cubic;

solver.xoff=1;
solver.yoff=1;
solver.zoff=1;

solver.temin=2.0;
solver.temax=12.0;
solver.sigmat=0.25;
solver.teshape=Te3dModel.XY_Gauss;

solver.nemin=1.0;
solver.nemax=1.0;
solver.sigman=0.25;
solver.neshape=Ne3dModel.Const3d;

solver.display=0;

