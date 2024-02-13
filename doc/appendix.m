function appendix(ind)
% function appendix(ind)

%
% $Id: appendix.m,v 1.4 2011/03/26 12:56:38 patrick Exp $
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

if nargin<1,
  ind = 1;
end

switch ind
  case 1,
    filename='appendix.dat';
    n=[21:2:419];
	case 2,
	  filename='appendix2.dat';
		n=[421:2:819];
	case 3,
	  filename='appendix3.dat';
		n=[821:2:1219];
	case 4,
	  filename='appendix4.dat';
		n=[1221:2:1619];
	case 5,
	  filename='appendix5.dat';
		n=[1621:2:2019];
	case 6,
	  filename='appendix6.dat';
		n=[2021:2:2419];
	case 7,
	  filename='appendix7.dat';
		n=[2421:2:2819];
	case 8,
	  filename='appendix8.dat';
		n=[2821:2:3219];
end

Nmax=10;
if rem(length(n),Nmax)~=0
	error('length(n) is not multiple of Nmax')
end
N=reshape(n,Nmax,fix(length(n)/Nmax));

fid=fopen(filename,'w');

for I=1:size(N,2),

n=N(:,I)';

npk=zeros(3,length(n));
j=1;
for i=n,
	[P,K]=calcnmudpack(i);
	npk(1,j)=i;
	npk(2,j)=P;
	npk(3,j)=K;
	j=j+1;
end

[tmp,ii]=sort(npk(3,:));
ii=flipud(ii(:));
npk(1:3,:)=npk(1:3,ii);

if I==1, 
	fprintf(fid,'\\begin{tabular}{l');
	for i=1:length(n), fprintf(fid,'|r'); end
	fprintf(fid,'}\n');
	fprintf(fid,'\\hline\\hline\n');
end
fprintf(fid,'$n$ (%d--%d)', min(npk(1,:)),max(npk(1,:)));
for i=1:length(n), fprintf(fid,' & %d', npk(1,i)); end
fprintf(fid,'\\\\\n');
fprintf(fid,'\\hline\n');
fprintf(fid,'$k$ $[p]$ ');
for i=1:length(n), fprintf(fid,'& %d [%d]', npk(3,i), npk(2,i)); end
fprintf(fid,'\\\\\n');
if I~=size(N,2),
	fprintf(fid,'\\hline\n');
else
	fprintf(fid,'\\hline\\hline\n');
	fprintf(fid,'\\end{tabular}\n'); 
end

end;

fclose(fid);

