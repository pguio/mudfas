function test7(dim)
% function test7(dim)

%
% $Id: test7.m,v 1.2 2011/03/26 12:56:39 patrick Exp $
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
	logfile='test7_2d.log';
	global nx ny
	global a kx ky X Y x y phit
elseif dim==3,
	logfile='test7_3d.log';
	global nx ny nz
	global a kx ky kz X Y Z x y z phit
end

if exist(logfile)==2, unix(['rm ' logfile]); end
eval(['diary ' logfile]);

global sd

SD=[1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e0];

if dim==2,
	fid=fopen('test7_2d.dat','w');
elseif dim==3,
	 fid=fopen('test7_3d.dat','w');
end

for I=1:length(SD),

	sd=SD(I);

	figure(1)
	t=cputime;
	if dim==2,
		[phi1,rho,solver]=mgfas2d('fun7');
		Te1=Te2d(solver,solver.ngrid);
	elseif dim==3,
		[phi1,rho,solver]=mgfas3d('fun7');
		Te1=Te3d(solver,solver.ngrid);
	end
	fprintf(1,'Time used=%.0f s\n', cputime-t);

	t=cputime;
	if dim==2,
		[PHI,RHO,solver]=fmg2d('fun7');
		phi5=PHI{solver.ngrid}(2:end-1,2:end-1);
		Te5=Te2d(solver,solver.ngrid);
	elseif dim==3
		[PHI,RHO,solver]=fmg3d('fun7');
		phi5=PHI{solver.ngrid}(2:end-1,2:end-1,2:end-1);
		Te5=Te3d(solver,solver.ngrid);
	end
	fprintf(1,'Time used=%.0f s\n', cputime-t);

if 0

	figure(2)
	argcaxis='auto';
	if  dim==2,
		plotphi(phi1,solver,argcaxis,[3 2 1],'mgfas')
		plotphi(phi5,solver,argcaxis,[3 2 2],'mudfas')
		colormap(gray)
		orient landscape; set(gcf,'PaperOrientation','portrait');
		if exist('exportfig'),
			exportfig(gcf,'test7_1_2d','format','eps');
		else
			print -deps test7_1_2d;
		end
	elseif dim==3,
		plotphi3d(phi1,solver,argcaxis,[3 2 1],'mgfas')
		plotphi3d(phi5,solver,argcaxis,[3 2 2],'mudfas')
		colormap(gray)
		orient landscape; set(gcf,'PaperOrientation','portrait');
		if exist('exportfig'),
			exportfig(gcf,'test7_1_3d','format','eps');
		else
			print -deps test7_1_3d;
		end
	end

end

	fprintf(fid,'%.0e & ',sd);

	figure(3)
	argcaxis=[-7 -3];
%	argcaxis='auto';
	if dim==2
		[d2phi,expphi,err]=plotd2d(phi1,rho,solver,argcaxis,[3 2 1],'mgfas');
	elseif dim==3,
		[d2phi,expphi,err]=plotd3d(phi1,rho,solver,argcaxis,[3 2 1],'mgfas');
	end
	fprintf(fid,'%.1e & %.1e &',mean(err(:)),std(err(:)));
	if dim==2,
		[d2phi,expphi,err]=plotd2d(phi5,rho,solver,argcaxis,[3 2 2],'mudfas');
	elseif dim==3,
		[d2phi,expphi,err]=plotd3d(phi5,rho,solver,argcaxis,[3 2 2],'mudfas');
	end
	fprintf(fid,'%.1e & %.1e \\\\\n',mean(err(:)),std(err(:)));


if 0
	colormap(gray)
	orient landscape; set(gcf,'PaperOrientation','portrait');
	if  dim==2
		if exist('exportfig'),
			exportfig(gcf,'test7_2_2d','format','eps');
		else
			print -deps test7_2_2d;
		end
	elseif dim==3
		if exist('exportfig'),
			exportfig(gcf,'test7_2_3d','format','eps');
		else
			print -deps test7_2_3d;
		end
	end
end

end

fclose(fid);

diary off
