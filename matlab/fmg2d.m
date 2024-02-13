function [U,F,solver]=fmg2d(f,varargin)
% function [U,F,solver]=fmg2d(f,varargin)

%
% $Id: fmg2d.m,v 1.5 2011/03/26 12:56:39 patrick Exp $
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

% Load default paramaters
solver=solver2d_def;
% Argument parsing 
solver=solver_par_parsing(solver,varargin{:});
% Load rhs and u(b.c.) 
[u,f,solver]=feval([f '_2d'],solver);
% Discretize pde at each grid level
solver=initpde2(solver);
% Display solver
view_solver(solver);
%
ngrid=solver.ngrid;
kcycle=solver.kcycle;
nx=solver.nx;
ny=solver.ny;
% set u, f in U{}, F{} and adjust right hand side
[U{ngrid},F{ngrid}]=swk2(u,f,solver);

if solver.iguess==0, % no initial guess at finest grid level!

% transfer down to all grid levels
	for k=ngrid:-1:2,
		[U{k-1},F{k-1}]=trsfc2(U{k},F{k},k-1,solver);
	end
% adjust right hand side at all grid levels in case
% rhs or specified b.c. in phi or gbdy changed
	for k=1:ngrid,
		[U{k},F{k}]=adjmd2(U{k},F{k},k,solver);
	end
%plot([U{ngrid}(3,2:end-1)', -F{ngrid}(2,:)']), pause
% execute one full multigrid cycle
	U{1}=smooth2d(U{1},F{1},1,solver);

	for k=2:ngrid
		U{k}=prolon2(U{k-1},U{k},k,solver);
	 	fprintf(1,'Entering mg2d with k=%d %s\n',k, repmat(' ',1,40));
	 	U=mg2d(U,F,k,solver);
	end
	if solver.display, 
		plot2dUF(U,F); 
	end;

else,
% transfer down to all grid levels
	for k=ngrid:-1:2,
		[U{k-1},F{k-1}]=trsfc2(U{k},F{k},k-1,solver);
	end
% adjust rhs at finest grid level only
	[U{ngrid},F{ngrid}]=adjmd2(U{ngrid},F{ngrid},ngrid,solver);
end;

% execute maxcy more multigrid k cycles from finest level
epsk=solver.tolmax*mean(abs(F{ngrid}(:)));
d=lop2d(U{ngrid},ngrid,solver)-F{ngrid};
tk=mean(abs(d(:)))-epsk;
if tk<0.0
	fprintf(1,'tk=%.4e mean |defect|=%.4e std |defect|=%.4e\n', ...
	tk,mean(abs(d(:))),std(abs(d(:))))
else
	for i=1:solver.maxcy,
		fprintf(1,'Entering iter#%3d of mg2d with k=%d\n',i,ngrid);
		U=mg2d(U,F,ngrid,solver);
		d=lop2d(U{ngrid},ngrid,solver)-F{ngrid};
		tk=mean(abs(d(:)))-epsk;
		if tk<0.0
			fprintf(1,'tk=%.4e mean |defect|=%.4e std |defect|=%.4e\n', ...
				tk,mean(abs(d(:))),std(abs(d(:))))
			break;
		end
	end
end

if solver.display, 
	plot2dUF(U,F); 
	pause;
end;

plotdefect2d(U{ngrid}(2:end-1,2:end-1),F{ngrid},solver)
