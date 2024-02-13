function y=nldF(phi,ne,te)
% function y=nldF(phi,ne,te)

%
% $Id: nldF.m,v 1.5 2011/03/26 12:56:39 patrick Exp $
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


y=-(ne./te).*exp(phi./te);


if 0
% From Mamun and Cairns, Stability of solitary waves in a magnetized
% non-thermal plasma, J. Plasma Physics, 56, 175-185, 1996
alpha=0.2;
% Eq. 2.2a
beta=4*alpha/(1+3*alpha);
% Eq. 2.1c
% y=-ne.*(1-beta*phi./te+beta*(phi./te).^2).*exp(phi./te);
% dyd/dphi
y=-ne.*(-beta./te+beta*2.0*(phi./te)./te).*exp(phi./te) - ...
	ne.*(1-beta*phi./te+beta*(phi./te).^2)./te.*exp(phi./te);
end
