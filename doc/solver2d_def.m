function solver=solver2d_def
% function solver=solver2d_def

%
% $Id: solver2d_def.m,v 1.3 2011/03/26 12:56:38 patrick Exp $
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

solver=struct('dim',[],'nx',[],'ny',[], ...
'nxa',[],'nxb',[],'nyc',[],'nyd',[], ...
'xa',[],'xb',[],'yc',[],'yd',[], ...
'ixp',[],'jyq',[],'iex',[],'jey',[],'ngrid',[], ...
'tolmax',[],'delta',[],'maxcy',[],'iguess',[], ...
'kcycle',[],'iprer',[],'ipost',[],'intpol',[], ...
'xoff',[],'yoff',[], ...
'temin',[],'temax',[],'sigmat',[],'teshape',[], ...
'display',[],...
'nk',[],'bk',[],'sk',[]);

solver.dim=2;

solver.nx=141;
solver.ny=349;

fid=fopen('nxny.dat','w');
fprintf(fid,'$%d\\times%d$~points ',solver.nx,solver.ny);
fclose(fid);

solver.nxa=bc.dirichlet;
solver.nxb=bc.dirichlet;
solver.nyc=bc.mixed;
solver.nyd=bc.mixed;

% Parameters for mixed boundary condition 
% For example at x=xa
% d(phi(xa, y))/dx + alfxa* phi(xa, y) = gbdxa(y)
solver.gbdxa=0;
solver.gbdxb=0;
solver.gbdyc=0;
solver.gbdyd=0;
%
solver.alfxa=0;
solver.alfxb=0;
solver.alfyc=0;
solver.alfyd=0;

solver.xa=0;
solver.xb=solver.nx-1;
solver.yc=0;
solver.yd=solver.ny-1;

fid=fopen('omega2d.dat','w');
fprintf(fid,'$\\Omega=[%.0f,%.0f]\\times[%.0f,%.0f]$', ...
	solver.xa,solver.xb,solver.yc,solver.yd);
fclose(fid);

solver.tolmax=5e-4;
solver.delta=5e-1;
solver.maxcy=5;
solver.iguess=0;

solver.kcycle=schedule.W;
solver.iprer=2;
solver.ipost=1;
solver.intpol=interp.Linear;

solver.xoff=1;
solver.yoff=1;

solver.temin=2.0;
solver.temax=2.0;
solver.sigmat=0.25;
solver.teshape=Te2dModel.X_Gauss;

solver.nemin=1.0;
solver.nemax=1.0;
solver.sigman=0.25;
solver.neshape=Te2dModel.Const2d;


solver.display=0;
