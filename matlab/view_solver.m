function view_solver(solver)
% function view_solver(solver)

%
% $Id: view_solver.m,v 1.4 2011/03/26 12:56:40 patrick Exp $
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


fprintf(1,'*************\n');
fprintf(1,'SOLVER SET UP\n');
fprintf(1,'*************\n');

fprintf(1,'NX = %3d ==> IXP = %3d IEX = %3d\n',solver.nx,solver.ixp,solver.iex);
fprintf(1,'NY = %3d ==> JYQ = %3d JEY = %3d\n',solver.ny,solver.jyq,solver.jey);
if solver.dim==3 & ...
	isfield(solver,'nz') & isfield(solver,'kzr') & isfield(solver,'kez'),
	fprintf(1,'NZ = %3d ==> KZR = %3d KEZ = %3d\n', ...
		solver.nz,solver.kzr,solver.kez);
end

fprintf(1,'==> NGRID = %2d\n', solver.ngrid);
for k=1:solver.ngrid,
	if solver.dim==2 & length([solver.nk{k}])==2,
		fprintf(1,'Grid# %2d ==> NX = %3d NY = %3d\n',k,solver.nk{k}(1:2));
	elseif solver.dim==3 & length(solver.nk{k})==3,
		fprintf(1,'Grid# %2d ==> NX = %3d NY = %3d NZ = %3d\n',k,solver.nk{k}(1:3));
	end
end

fprintf(1,'XA = %5.1f ==> %s\n',solver.xa,upper(msgboundary(solver.nxa)));
fprintf(1,'XB = %5.1f ==> %s\n',solver.xb,upper(msgboundary(solver.nxb)));
fprintf(1,'YC = %5.1f ==> %s\n',solver.yc,upper(msgboundary(solver.nyc)));
fprintf(1,'YD = %5.1f ==> %s\n',solver.yd,upper(msgboundary(solver.nyd)));
if solver.dim==3 & isfield(solver,'nze'),
	fprintf(1,'ZE = %5.1f ==> %s\n',solver.ze,upper(msgboundary(solver.nze)));
end
if solver.dim==3 & isfield(solver,'nzf'),
	fprintf(1,'ZF = %5.1f ==> %s\n',solver.zf,upper(msgboundary(solver.nzf)));
end

fprintf(1,'TOLMAX = %.2E\n',solver.tolmax);
fprintf(1,'DELTA  = %.2E\n',solver.delta);
fprintf(1,'MAXCY  = %3d\n',solver.maxcy);
fprintf(1,'IGUESS = %3d\n',solver.iguess);
fprintf(1,'KCYCLE = %3d ==> %s\n',solver.kcycle,upper(msgcycle(solver.kcycle)));
fprintf(1,'IPRER  = %3d\n',solver.iprer);
fprintf(1,'IPOST  = %3d\n',solver.ipost);
fprintf(1,'INTPOL = %3d ==> %s\n',solver.intpol,upper(msgintpol(solver.intpol)));

fprintf(1,'TEMIN   = %4.2f\n',solver.temin);
fprintf(1,'TEMAX   = %4.2f\n',solver.temax);
fprintf(1,'SIGMAT  = %4.2f\n',solver.sigmat);
if solver.dim==2 
	fprintf(1,'TESHAPE = %s\n', msgshape2d(solver.teshape));
elseif solver.dim==3
	fprintf(1,'TESHAPE = %s\n', msgshape3d(solver.teshape));
end

fprintf(1,'NEMIN   = %4.2f\n',solver.nemin);
fprintf(1,'NEMAX   = %4.2f\n',solver.nemax);
fprintf(1,'SIGMAN  = %4.2f\n',solver.sigman);
if solver.dim==2 
	fprintf(1,'NESHAPE = %s\n', msgshape2d(solver.neshape));
elseif solver.dim==3
	fprintf(1,'NESHAPE = %s\n', msgshape3d(solver.neshape));
end

fprintf(1,'*************\n');


function msg=msgboundary(BC)
solver_const;

if (BC~=bc.periodic & BC~=bc.dirichlet & BC~=bc.mixed)
	error('Not a valid boundary condition.\n');
end
msg=BCTypeName{BC+1};


function msg=msgcycle(kcycle)
solver_const;

if kcycle<=0
	error('Not a valid schedule, kcycle must be >0.');
elseif kcycle~=schedule.V & kcycle~=schedule.W
	msg='UNKNOWN PREDEFINED CYCLE';
else
	msg=MGScheduleName{kcycle+1};
end


function msg=msgintpol(intpol)
solver_const;

if intpol~=interp.Linear & intpol~=interp.Cubic
	error('Not a valid interpolation operator.');
end
msg=InterpTypeName{intpol};

function msg=msgshape2d(shape)
solver_const;


if shape<0 | shape>length(struct2cell(Te2dModel))
	error('Not a valid shape.')
end
msg=Shape2DName{shape};

function msg=msgshape3d(shape)
solver_const;

if shape<0 | shape>length(struct2cell(Te3dModel))
	error('Not a valid shape.')
end
msg=Shape3DName{shape};

