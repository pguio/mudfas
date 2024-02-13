function solver=solver_par_parsing(solver,varargin)
% function solver=solver_par_parsing(solver,varargin)

%
% $Id: solver_par_parsing.m,v 1.3 2011/03/26 12:56:39 patrick Exp $
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

if isempty(varargin),
	return
end

argc=length(varargin);

if rem(argc,2)~=0,
	error('Extra arguments should be in pair parameter/value')
end

for i=1:2:argc,

	param=varargin{i};
	val=varargin{i+1};
	switch lower(param)
		case 'nx', solver.nx=val;
		case 'ny', solver.ny=val;
		case 'nz', solver.nz=val;

		case 'nxa', solver.nxa=val;
		case 'nxb', solver.nxb=val;
		case 'nyc', solver.nyc=val;
		case 'nyd', solver.nyd=val;
		case 'nze', solver.nze=val;
		case 'nzf', solver.nzf=val;

		case 'xa', solver.xa=val;
		case 'xb', solver.xb=val;
		case 'yc', solver.yc=val;
		case 'yd', solver.yd=val;
		case 'ze', solver.ze=val;
		case 'zf', solver.zf=val;

		case 'gbdxa', solver.gbdxa=val;
		case 'gbdxb', solver.gbdxb=val;
		case 'gbdyc', solver.gbdyc=val;
		case 'gbdyd', solver.gbdyd=val;
		case 'gbdze', solver.gbdze=val;
		case 'gbdzf', solver.gbdzf=val;

		case 'alfxa', solver.alfxa=val;
		case 'alfxb', solver.alfxb=val;
		case 'alfyc', solver.alfyc=val;
		case 'alfyd', solver.alfyd=val;
		case 'alfze', solver.alfze=val;
		case 'alfzf', solver.alfzf=val;

		case 'tolmax', solver.tolmax=val;
		case 'delta', solver.delta=val;
		case 'maxcy', solver.maxcy=val;
		case 'iguess', solver.iguess=val;

		case 'kcycle', solver.kcycle=val;
		case 'iprer', solver.iprer=val;
		case 'ipost', solver.ipost=val;
		case 'intpol', solver.intpol=val;

		case 'temin', solver.temin=val;
		case 'temax', solver.temax=val;
		case 'sigmat', solver.sigmat=val;
		case 'teshape', solver.teshape=val;

		case 'nemin', solver.nemin=val;
		case 'nemax', solver.nemax=val;
		case 'sigman', solver.sigman=val;
		case 'neshape', solver.neshape=val;

		case 'display', solver.display=val;

	end

end

