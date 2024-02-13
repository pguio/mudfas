function test3(dim)
% function test3(dim)

%
% $Id: test3.m,v 1.2 2011/03/26 12:56:39 patrick Exp $
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
	logfile='test3_2d.log';
elseif dim==3,
	logfile='test3_3d.log';
end

if exist(logfile)==2, unix(['rm ' logfile]); end
eval(['diary ' logfile]);

KY=[2:2:12]*pi;

if dim==2,
	global a kx ky X Y x y phit
elseif dim==3
	global a kx ky kz X Y Z x y z phit
end

if dim==2,
	fid=fopen('test3_2d.dat','w');
elseif dim==3,
	fid=fopen('test3_3d.dat','w');
end

for I=1:6,

ky=KY(I);

figure(1)
t=cputime;
if dim==2,
	[phi1,rho,solver]=mgfas2d('fun3');
	Te1=Te2d(solver,solver.ngrid);
elseif dim==3
	[phi1,rho,solver]=mgfas3d('fun3');
	Te1=Te3d(solver,solver.ngrid);
end
fprintf(1,'Time used=%.0f s\n', cputime-t);

t=cputime;
if dim==2,
	[PHI,RHO,solver]=fmg2d('fun3');
	phi5=PHI{solver.ngrid}(2:end-1,2:end-1);
	Te5=Te2d(solver,solver.ngrid);
elseif dim==3
	[PHI,RHO,solver]=fmg3d('fun3');
	phi5=PHI{solver.ngrid}(2:end-1,2:end-1,2:end-1);
	Te5=Te3d(solver,solver.ngrid);
end
fprintf(1,'Time used=%.0f s\n', cputime-t);

xa=solver.xa; xb=solver.xb;
yc=solver.yc; yd=solver.yd;
if dim==3,
  ze=solver.ze; zf=solver.zf;
end

if dim==2,
	phit = a./(kx^2/(xb-xa)^2+ky^2/(yd-yc)^2+1./Te1).*phit;
elseif dim==3,
	phit = a./(kx^2/(xb-xa)^2+ky^2/(yd-yc)^2+kz^2/(zf-ze)^2+1./Te1).*phit;
end

if 0

figure(2)
argcaxis='auto';
subplot(3,3,I);
ys=(yd-yc)*y;
ii=fix(solver.nx/2);
if dim==2,
	plot(ys,phit(ii,:),'-',ys,phi5(ii,:),'--',ys,phi1(ii,:),'-.');
elseif dim==3,
	kk=fix(solver.nz/2);
	plot(ys,phit(ii,:,kk),ys,phi5(ii,:,kk),'--',ys,phi1(ii,:,kk),'-.');
end
set(gca,'xlim',[min(ys) max(ys)]);
set(gca,'ylim',[min([phit(:);phi1(:);phi5(:)]) max([phit(:);phi1(:);phi5(:)])]);
title(sprintf('\\phi(y) with ky=%.0f\\pi',ky/pi))
orient landscape; set(gcf,'PaperOrientation','portrait');
if dim==2,
	if exist('exportfig'),
		exportfig(gcf,'test3_1_2d','format','eps');
	else
		print -deps test3_1_2d;
	end
elseif dim==3,
	if exist('exportfig'),
		exportfig(gcf,'test3_1_3d','format','eps');
	else
		print -deps test3_1_3d;
	end
end

end

fprintf(fid,'$%.0f\\pi$ & ',ky/pi);

figure(3)
argcaxis='auto';
if dim==2,
	[d2phi,expphi,err]=plotd2d(phit,rho,solver,argcaxis,[3 3 1],'lin');
elseif dim==3,
	[d2phi,expphi,err]=plotd3d(phit,rho,solver,argcaxis,[3 3 1],'lin');
end
fprintf(fid,'%.1e & %.1e & ',mean(err(:)),std(err(:)));
if dim==2,
	[d2phi,expphi,err]=plotd2d(phi1,rho,solver,argcaxis,[3 3 2],'mgfas');
elseif dim==3,
	[d2phi,expphi,err]=plotd3d(phi1,rho,solver,argcaxis,[3 3 2],'mgfas');
end
fprintf(fid,'%.1e & %.1e & ',mean(err(:)),std(err(:)));
if dim==2,
	[d2phi,expphi,err]=plotd2d(phi5,rho,solver,argcaxis,[3 3 3],'mudfas');
elseif dim==3,
	[d2phi,expphi,err]=plotd3d(phi5,rho,solver,argcaxis,[3 3 3],'mudfas');
end
fprintf(fid,'%.1e & %.1e \\\\\n',mean(err(:)),std(err(:)));

end

fclose(fid);

diary off
