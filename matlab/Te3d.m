function te=Te3d(solver,k)
% function te=Te3d(solver,k)

%
% $Id: Te3d.m,v 1.5 2011/03/26 12:56:39 patrick Exp $
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
nz=solver.nk{k}(3);


switch solver.teshape
	case Te3dModel.XY_Gauss,
		[X,Y,Z]=ndgrid(linspace(-1,1,nx),linspace(-1,1,ny),zeros(1,nz));
		te=(temax-temin)*exp(-(X.^2+Y.^2)/2/sigmat^2)+temin;
		return
	case Te3dModel.XZ_Gauss,
		[X,Y,Z]=ndgrid(linspace(-1,1,nx),zeros(1,ny),linspace(-1,1,nz));
		te=(temax-temin)*exp(-(X.^2+Z.^2)/2/sigmat^2)+temin;
		return
	case Te3dModel.Const3d,
		te=repmat(temin,[nx,ny,nz]);
		return
	otherwise
		error('Not a valid Te shape');
end		


