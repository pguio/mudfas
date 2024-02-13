function [u,f,solver]=fun8(solver)
% function [u,f,solver]=fun8(solver)

%
% $Id: fun8_2d.m,v 1.2 2011/03/26 12:56:39 patrick Exp $
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

u=zeros(solver.nx,solver.ny);

x=linspace(0,solver.nx-1,solver.nx);
y=linspace(0,solver.ny-1,solver.ny);
[X,Y]=ndgrid(x,y);

K=1e-2;
L=solver.ny-1;
n0=K*L/(1-exp(-K*L));
Te=1;
f=n0*exp(-K*Y);

fprintf('n0 range=(%f,%f)\n', n0, n0*exp(-K*L));

solver.temin=Te;
solver.temax=solver.temin;


y0=1/K*log(n0);

solver.analytic_u=-K*Te*(Y-y0);

if 0,
	solver.gbdyc=-K*Te;
	solver.gbdyd=-K*Te;
else
	u([1, end],:)=solver.analytic_u([1, end],:);
	u(:,[1, end])=solver.analytic_u(:,[1, end]);
	solver.nxa=1;
	solver.nxb=1;
	solver.nyc=1;
	solver.nyd=1;
end
