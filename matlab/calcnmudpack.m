function [p,k]=calcnmudpack(n)
% function [p,k]=calcnmudpack(n)

%
% $Id: calcnmudpack.m,v 1.3 2011/03/26 12:56:39 patrick Exp $
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

if 1

P=[2 3:2:fix(n/2)];
K=1:fix(n/2)+1;

[P,K]=ndgrid(P,K);
N=bitshift(P,K-1)+1;
[i,j]=find(N==n);
if ~isempty(i),
	p=P(i,j);
	k=K(i,j);
else
	error('calnmudpack: p,k not found');
end
	

else

%for p=[primes(fix(n/2)) fix(n/2)],
for p=[2 3:2:fix(n/2)],
	for k=1:fix(n/2)+1,
		%N=fix(p*2^(k-1)+1);
		N=bitshift(p,k-1)+1;
		%fprintf(1,'p=%2.d, k=%2.d, n=%3.d\n', p, k, N);
		if (N>n),
		  break;
		end;
		if (N==n)
%		  fprintf(1,'p=%3.d, k=%3.d, n=%3.d\n', p, k, N);
		  return; 
		end;
	end
end

error('calnmudpack: p,k not found');

end
