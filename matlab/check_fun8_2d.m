function check_fun8_2d
% function check_fun8_2d

%
% $Id: check_fun8_2d.m,v 1.3 2011/03/26 12:56:39 patrick Exp $
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

[u,f,s]=fmg2d('fun8');

U=mean(u{3}(2:end-1,2:end-1),1);
analytic_u=mean(s.analytic_u,1);
ni=mean(-f{3},1);
ne=ones(size(ni));

x=linspace(0,s.ny-1,s.ny);

subplot(211)
plot(x,ni,'-',x,ne,'--')
title('density');
legend('ion','electron');

subplot(212);
plot(x,U,'-',x,analytic_u,'--')
title('electric potential')
legend('numerical','analytical');
