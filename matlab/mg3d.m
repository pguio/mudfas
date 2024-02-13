function U=mg3d(U,F,k,solver)
% function U=mg3d(U,F,k,solver)

%
% $Id: mg3d.m,v 1.5 2011/03/26 12:56:39 patrick Exp $
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

if solver.display>=0,
	s=sprintf('%s>%s In  mg3d with k=%2d', ...
		repmat('-',1,2*k), repmat(' ',1,2*(solver.ngrid-k+1)), k);
	fprintf(1,'%s%s', s, char(repmat(8,1,length(s))));
end

iprer=solver.iprer;
ipost=solver.ipost;
kcycle=solver.kcycle;

if k==1,
  U{1}=smooth3d(U{1},F{1},1,solver);
else,
% pre-smoothing
	for i=1:iprer,
  	U{k}=smooth3d(U{k},F{k},k,solver);
	end
% form residual=-defect
	rk=F{k}-lop3d(U{k},k,solver);
% choose u{k-1}
	U{k-1}=resu3(U{k},U{k-1},k-1,solver);
% 
	F{k-1}=lop3d(U{k-1},k-1,solver)+solver.sk{k-1}* ...
		res3(rk,zeros(size(F{k-1})),k-1,solver);

% adaptive criteria
	tkm1=mean(abs(rk(:)))-solver.sk{k}*mean(abs(F{k}(:)));
	epskm1=solver.delta*solver.sk{k-1}*mean(abs(rk(:)));

	Uoldkm1=U{k-1};

	for i=1:kcycle,
		U=mg3d(U,F,k-1,solver);
% adaptive criteria
		t=lop3d(U{k-1},k-1,solver)-F{k-1};
		tkm1=mean(abs(t(:)))-epskm1;
		if tkm1<0, 
			break; 
		end;
	end

% inject correction
	U{k}=U{k}+1/solver.sk{k-1}* ...
		prolon3(U{k-1}-Uoldkm1,zeros(size(U{k})),k,solver);

% post-smoothing
	for i=1:ipost,
  	U{k}=smooth3d(U{k},F{k},k,solver);
	end
end

if solver.display>=0,
	s=sprintf('<%s%s Out mg3d with k=%2d', ...
		repmat('-',1,2*k),repmat(' ',1,2*(solver.ngrid-k+1)),k);
	fprintf(1,'%s%s', s, char(repmat(8,1,length(s))));
end
