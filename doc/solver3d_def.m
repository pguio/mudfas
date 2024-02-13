function solver=solver3d_def
% function solver=solver3d_def

%
% $Id: solver3d_def.m,v 1.3 2011/03/26 12:56:39 patrick Exp $
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

solver.nx=45;
solver.ny=85;
solver.nz=45;

fid=fopen('nxnynz.dat','w');
fprintf(fid,'$%d\\times%d\\times%d$~points ',solver.nx,solver.ny,solver.nz);
fclose(fid);

solver.nxa=bc.dirichlet;
solver.nxb=bc.dirichlet;
solver.nyc=bc.mixed;
solver.nyd=bc.mixed;
solver.nze=bc.dirichlet;
solver.nzf=bc.dirichlet;

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

solver.xa=0;
solver.xb=solver.nx-1;
solver.yc=0;
solver.yd=solver.ny-1;
solver.ze=0;
solver.zf=solver.nz-1;

fid=fopen('omega3d.dat','w');
fprintf(fid,'$\\Omega=[%.0f,%.0f]\\times[%.0f,%.0f]\\times[%.0f,%.0f]$.', ...
  solver.xa,solver.xb,solver.yc,solver.yd,solver.ze,solver.zf);
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
solver.zoff=1;

solver.temin=2.0;
solver.temax=2.0;
solver.sigmat=0.25;
solver.teshape=Te3dModel.XZ_Gauss;

solver.nemin=1.0;
solver.nemax=1.0;
solver.sigman=0.25;
solver.neshape=Te3dModel.Const3d;

solver.display=0;

