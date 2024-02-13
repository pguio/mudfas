function te=Te2d(solver,k)
% function te=Te2d(solver,k)

%
% $Id: Te2d.m,v 1.5 2011/03/26 12:56:39 patrick Exp $
%
% Copyright (c) 2000-2011 Patrick Guio <patrick.guio@gmail.com>
% All Rights Reserved.
%
% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the Free
% Software Foundation; either version 2.  of the License, or (at your
% option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
% or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
% for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <http://www.gnu.org/licenses/>.


solver_const;

temin=solver.temin;
temax=solver.temax;
sigmat=solver.sigmat;
nx=solver.nk{k}(1);
ny=solver.nk{k}(2);

switch solver.teshape
	case Te2dModel.X_Gauss
		[X,Y]=ndgrid(linspace(-1,1,nx),zeros(1,ny));
		te=(temax-temin)*exp(-X.^2/2/sigmat^2)+temin;
		return
	case Te2dModel.Y_Gauss
		[X,Y]=ndgrid(zeros(1,nx),linspace(-1,1,ny));
		te=(temax-temin)*exp(-Y.^2/2/sigmat^2)+temin;
		return
	case Te2dModel.X_Dep
		[X,Y]=ndgrid(zeros(1,nx),linspace(-1,1,ny));
		te=(temax-temin)*(1.0-exp(-Y.^2/2/sigmat^2))+temin;
		return
	case Te2dModel.X_Erf
		[X,Y]=ndgrid(zeros(1,nx),linspace(-1,1,ny));
		te=0.5*(temax-temin)*(1.0+erf(Y/sigmat))+temin;
		return
	case Te2dModel.X_Atan
		[X,Y]=ndgrid(zeros(1,nx),linspace(-1,1,ny));
		te=0.5*(temax-temin)*(2.0*atan(Y/sigmat)/pi+1.0)+temin;
		return
	case Te2dModel.Const2d
		te=repmat(temin,[nx,ny]);
		return
	otherwise
		error('Not a valid Te shape');
end		

