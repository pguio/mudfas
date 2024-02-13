function calctable(n,filename)
% function calctable(n,filename)

%
% $Id: calctable.m,v 1.2 2011/03/26 12:56:38 patrick Exp $
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

n=n+1;
npk=zeros(3,length(n));
j=1;
for i=n,
	[P,K]=calcnmudpack(i);
	npk(1,j)=i;
	npk(2,j)=P;
	npk(3,j)=K;
	j=j+1;
end

fid=fopen(filename,'w');

fprintf(fid,'\\begin{tabular}{r');
for i=1:length(n), fprintf(fid,'|r'); end
fprintf(fid,'}\n');
fprintf(fid,'\\hline\\hline\n');
fprintf(fid,'$n$');
for i=1:length(n), fprintf(fid,' & %d', npk(1,i)); end
fprintf(fid,'\\\\\n');
fprintf(fid,'\\hline\n');
%fprintf(fid,'$k$ $(p)$ ');
fprintf(fid,'$k$ ');
%for i=1:length(n), fprintf(fid,'& %d (%d)', npk(3,i),  npk(2,i)); end
for i=1:length(n), fprintf(fid,'& %d ', npk(3,i)); end
fprintf(fid,'\\\\\n');
fprintf(fid,'\\hline\\hline\n');
fprintf(fid,'\\end{tabular}');

fclose(fid);

