function plot2dUF(U,F)
% function plot2dUF(U,F)

%
% $Id: plot2dUF.m,v 1.5 2011/03/26 12:56:39 patrick Exp $
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

	subplot(2,nb,i),
	if 1
		imagesc(U{i}(2:end-1,2:end-1)'); axis xy
		colorbar('h');
	else
		pcolor(U{i}(2:end-1,2:end-1)); 
		shading flat
	end
	title('\phi');
	drawnow

	subplot(2,nb,i+nb),
	if 1
		imagesc(F{i}'); axis xy
		colorbar('h');
	else
		pcolor(F{i});
		shading flat
	end
	title('F');
	drawnow

end
