function test6(dim)
% function test6(dim)

%
% $Id: test6.m,v 1.2 2011/03/26 12:56:39 patrick Exp $
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
	logfile='test6_2d.log';
	global nx ny
	global a kx ky X Y x y phit
elseif dim==3,
	logfile='test6_3d.log';
	global nx ny nz
	global a kx ky kz X Y Z x y z phit
end

if exist(logfile)==2, unix(['rm ' logfile]); end
eval(['diary ' logfile]);

if dim==2,
	NX=[143 141 137 145 161 193];
	NY=[351 349 345 337 353 321];
	calctable(NX-1,'npkx2d.dat');
	calctable(NY-1,'npky2d.dat');
elseif dim==3,
	NX=[43 45 41 49];
	NY=[83 85 89 81];
	NZ=[43 45 41 49];
	calctable(NX-1,'npkx3d.dat');
	calctable(NY-1,'npky3d.dat');
end

if dim==2,
	fid=fopen('test6_2d.dat','w');
elseif dim==3,
	 fid=fopen('test6_3d.dat','w');
end

for I=1:length(NX),

	nx=NX(I);
	ny=NY(I);
	if dim==3,
		nz=NZ(I);
	end

	figure(1)
	t=cputime;
	if dim==2,
		[phi1,rho,solver]=mgfas2d('fun6');
		Te1=Te2d(solver,solver.ngrid);
	elseif dim==3,
		[phi1,rho,solver]=mgfas3d('fun6');
		Te1=Te3d(solver,solver.ngrid);
	end
	fprintf(1,'Time used=%.0f s\n', cputime-t);

	t=cputime;
	if dim==2,
		[PHI,RHO,solver]=fmg2d('fun6');
		phi5=PHI{solver.ngrid}(2:end-1,2:end-1);
		Te5=Te2d(solver,solver.ngrid);
	elseif dim==3
		[PHI,RHO,solver]=fmg3d('fun6');
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
			exportfig(gcf,'test6_1_2d','format','eps');
		else
			print -deps test6_1_2d;
		end
	elseif dim==3,
		plotphi3d(phi1,solver,argcaxis,[3 2 1],'mgfas')
		plotphi3d(phi5,solver,argcaxis,[3 2 2],'mudfas')
		colormap(gray)
		orient landscape; set(gcf,'PaperOrientation','portrait');
		if exist('exportfig'),
			exportfig(gcf,'test6_1_3d','format','eps');
		else
			print -deps test6_1_3d;
		end
	end

	end

	if  dim==2,
		fprintf(fid,'$%d\\times%d$ $(%d)$ & ',nx,ny,length(PHI));
	elseif dim==3,
		fprintf(fid,'$%d\\times%d\\times%d$  $(%d)$ & ',nx,ny,nz,length(PHI));
	end

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

	colormap(gray)
	orient landscape; set(gcf,'PaperOrientation','portrait');
	if  dim==2
		if exist('exportfig'),
			exportfig(gcf,'test6_2_2d','format','eps');
		else
			print -deps test6_2_2d;
		end
	elseif dim==3
		if exist('exportfig'),
			exportfig(gcf,'test6_2_3d','format','eps');
		else
			print -deps test6_2_3d;
		end
	end

end

fclose(fid);

diary off
