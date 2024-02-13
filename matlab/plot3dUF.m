function plot3dUF(U,F)
% function plot3dUF(U,F)

%
% $Id: plot3dUF.m,v 1.5 2011/03/26 12:56:39 patrick Exp $
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

nb=length(U);

nx=fix(sqrt(nb));
ny=fix(nb/nx);

for i=1:nb,

	u=U{i}(2:end-1,2:end-1,2:end-1);
	if 1
		subplot(6,nb,i),
		imagesc(squeeze(u(:,:,fix(end/2)))'), axis xy, title('\phi XY')
		subplot(6,nb,i+nb),
		imagesc(squeeze(u(:,fix(end/2),:))'), axis xy, title('\phi XZ')
		subplot(6,nb,i+2*nb),
		imagesc(squeeze(u(fix(end/2),:,:))'), axis xy, title('\phi YZ')
	else
		nx=size(u,1);
		ny=size(u,2);
		nz=size(u,3);
		if min(U{i}(:))~=max(U{i}(:)),
			subplot(2,nb,i),
			fprintf(1,'sx=%f sy=%f sz=%f\n',nx/2,ny/2,nz/2);
			slice(u,nx/2,ny/2,nz/2);
			shading flat;
			set(gca,'xlim',[1 nx],'ylim',[1 ny],'zlim',[1 nz]);
		end
	end
	drawnow

	if 1
		subplot(6,nb,i+3*nb),
		imagesc(squeeze(F{i}(:,:,fix(end/2)))'), axis xy, title('F XY')
		subplot(6,nb,i+4*nb),
		imagesc(squeeze(F{i}(:,fix(end/2),:))'), axis xy, title('F XZ')
		subplot(6,nb,i+5*nb),
		imagesc(squeeze(F{i}(fix(end/2),:,:))'), axis xy, title('F YZ')
	else
		nx=size(F{i},1);
		ny=size(F{i},2);
		nz=size(F{i},3);
		if min(F{i}(:))~=max(F{i}(:)),
			subplot(2,nb,i+nb),
			fprintf(1,'sx=%f sy=%f sz=%f\n',nx/2,ny/2,nz/2);
			slice(F{i},nx/2,ny/2,nz/2);
			set(gca,'xlim',[1 nx],'ylim',[1 ny],'zlim',[1 nz]);
			shading flat;
		end
	end
	drawnow

end

