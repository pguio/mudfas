function [u,f,solver]=fun4(solver)
% function [u,f,solver]=fun4(solver)

%
% $Id: fun4_3d.m,v 1.3 2011/03/26 12:56:39 patrick Exp $
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

sd_id=hdfsd('start','beamBS3d.hdf','rdonly');

attr_idx = hdfsd('findattr',sd_id,'lx');
[lx,status] = hdfsd('readattr',sd_id,attr_idx);

attr_idx = hdfsd('findattr',sd_id,'ly');
[ly,status] = hdfsd('readattr',sd_id,attr_idx);

attr_idx = hdfsd('findattr',sd_id,'lz');
[lz,status] = hdfsd('readattr',sd_id,attr_idx);

idx = hdfsd('nametoindex',sd_id','rho');
sds_id = hdfsd('select',sd_id,idx);
[name,rank,dimsizes,data_type,nattrs,status] = hdfsd('getinfo',sds_id);

[rho,status]=hdfsd('readdata',sds_id,[0 0 0],[],[dimsizes(1:3)]);


hdfsd('endaccess',sds_id);
hdfsd('end',sd_id);

solver_const;
if 1
	rho=shiftdim(rho,2);
	solver.xa=0;
	solver.xb=lz;
	solver.yc=0;
	solver.yd=lx;
	solver.ze=0;
	solver.zf=ly;
	solver.teshape=Te3dModel.XY_Gauss;
else
	solver.xa=0;
	solver.xb=lx;
	solver.yc=0;
	solver.yd=ly;
	solver.ze=0;
	solver.zf=lz;
	solver.teshape=Te3dModelXZ_Gauss;
end

solver.nx=size(rho,1);
solver.ny=size(rho,2);
solver.nz=size(rho,3);

%solver.temax=solver.temin;

u=zeros(solver.nx,solver.ny,solver.nz);

f=double(rho);


