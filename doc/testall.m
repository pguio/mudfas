function testall
% function testall

%
% $Id: testall.m,v 1.2 2011/03/26 12:56:39 patrick Exp $
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

for i=[1 3:8],
	disp(['test' int2str(i) '(2)']);
	eval(['test' int2str(i) '(2)']);
end

for i=[1 3:8],
	disp(['test' int2str(i) '(3)']);
	eval(['test' int2str(i) '(3)']);
end

close all

appendix


