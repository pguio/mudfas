function s = solver_const
% s = solver_const

%
% $Id: solver_const.m,v 1.6 2011/03/26 12:56:39 patrick Exp $
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

% Boundary conditions stuff
s.bc=struct('periodic',0,'dirichlet',1,'mixed',2);

s.BCTypeName={...
'PERIODIC BOUNDARY',...
'VALUE SPECIFIED AT BOUNDARY (DIRICHLET)',...
'MIXED DERIVATIVE BOUNDARY   (DIRICHLET/NEUMANN)'};

% Multi-grid schedule stuff
s.schedule=struct('V',1,'W',2);

s.MGScheduleName={...
'','V CYCLE','W CYCLE'};

% Interpolation operator stuff
s.interp=struct('Linear',1,'Cubic',3);

s.InterpTypeName={...
'LINEAR INTERPOLATION','','CUBIC INTERPOLATION'};

% 2D electron temperature model
s.Te2dModel=struct('Y_Gauss',1,'X_Gauss',2,'X_Dep',3,'X_Erf',4,...
	'X_Atan',5,'Const2d',6);

% 2D electron density model
s.Ne2dModel=struct('Y_Gauss',1,'X_Gauss',2,'X_Dep',3,'X_Erf',4,...
  'X_Atan',5,'Const2d',6);

s.Shape2DName={...,
'Y-GAUSSIAN',...
'X-GAUSSIAN',...
'X-DEPLETION',...
'X-ERF',...
'X-ATAN',...
'CONST'};

% 3D electron temperature model
s.Te3dModel=struct('XY_Gauss',1,'XZ_Gauss',2,'Const3d',3);

% 3D electron temperature model
s.Ne3dModel=struct('XY_Gauss',1,'XZ_Gauss',2,'Const3d',3);

s.Shape3DName={...,
'XY-GAUSSIAN', ...
'XZ-GAUSSIAN', ...
'CONST'};




