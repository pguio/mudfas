function [u,f,solver]=fun2(solver)
% function [u,f,solver]=fun2(solver)

%
% $Id: fun10_2d.m,v 1.4 2011/03/26 12:56:39 patrick Exp $
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

solver.nx=257;
solver.ny=257;

%solver.nx=505;
%solver.ny=505;

solver.xa=-60;
solver.xb=60;
solver.yc=-40;
solver.yd=80;


u=zeros(solver.nx,solver.ny);

f=ones(solver.nx,solver.ny);

[x,y]=meshgrid(...
	linspace(solver.yc,solver.yd,solver.ny), ...
	linspace(solver.xa,solver.xb,solver.nx)  ...
	);

R=5;
%phiOverTe=-4;
%phiOverTe=-2;
phiOverTe=-1.4;

ii=find(sqrt(x.^2+y.^2)<=R);
%ii=find(abs(x)<=R & abs(y)<=R);


f(ii)=exp(phiOverTe);

