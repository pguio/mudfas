function test8(dim)
% function test8(dim)

%
% $Id: test8.m,v 1.2 2011/03/26 12:56:39 patrick Exp $
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

close all

if dim==2,
	logfile='test8_2d.log';
	global nx ny
	global a kx ky X Y x y phit
elseif dim==3,
	logfile='test8_3d.log';
	global nx ny nz
	global a kx ky kz X Y Z x y z phit
end

if exist(logfile)==2, unix(['rm ' logfile]); end
eval(['diary ' logfile]);


global iprer ipost

IPRER=[1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4];
IPOST=[1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4];

if dim==2,
	fid=fopen('test8_2d.dat','w');
elseif dim==3,
	 fid=fopen('test8_3d.dat','w');
end

for I=1:length(IPRER),

	iprer=IPRER(I);
	ipost=IPOST(I);

	figure(1)
	t=cputime;
	if dim==2,
		[PHI,RHO,solver]=fmg2d('fun8');
		phi5=PHI{solver.ngrid}(2:end-1,2:end-1);
		rho=RHO{solver.ngrid};
		Te5=Te2d(solver,solver.ngrid);
	elseif dim==3
		[PHI,RHO,solver]=fmg3d('fun8');
		phi5=PHI{solver.ngrid}(2:end-1,2:end-1,2:end-1);
		rho=RHO{solver.ngrid};
		Te5=Te3d(solver,solver.ngrid);
	end
	fprintf(1,'Time used=%.0f s\n', cputime-t);

if 0

	figure(2)
	argcaxis='auto';
	if  dim==2,
		plotphi(phi5,solver,argcaxis,[3 1 1],'mudfas')
		colormap(gray)
		orient landscape; set(gcf,'PaperOrientation','portrait');
		if exist('exportfig'),
			exportfig(gcf,'test8_1_2d','format','eps');
		else
			print -deps test8_1_2d;
		end
	elseif dim==3,
		plotphi3d(phi5,solver,argcaxis,[3 1 1],'mudfas')
		colormap(gray)
		orient landscape; set(gcf,'PaperOrientation','portrait');
		if exist('exportfig'),
			exportfig(gcf,'test8_1_3d','format','eps');
		else
			print -deps test8_1_3d;
		end
	end

end

	fprintf(fid,'(%d, %d) & ', iprer, ipost);

	figure(3)
	argcaxis=[-7 -3];
%	argcaxis='auto';
	if dim==2,
		[d2phi,expphi,err]=plotd2d(phi5,rho,solver,argcaxis,[3 1 1],'mudfas');
	elseif dim==3,
		[d2phi,expphi,err]=plotd3d(phi5,rho,solver,argcaxis,[3 1 1],'mudfas');
	end
	fprintf(fid,'%.1e & %.1e \\\\\n',mean(err(:)),std(err(:)));


if 0
	colormap(gray)
	orient landscape; set(gcf,'PaperOrientation','portrait');
	if  dim==2
		if exist('exportfig'),
			exportfig(gcf,'test8_2_2d','format','eps');
		else
			print -deps test8_2_2d;
		end
	elseif dim==3
		if exist('exportfig'),
			exportfig(gcf,'test8_2_3d','format','eps');
		else
			print -deps test8_2_3d;
		end
	end
end

end

fclose(fid);

diary off
