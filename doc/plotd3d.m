function [d2phi,expphi,err]=plotd3d(u,f,solver,argcaxis,mnp,func)
% function [d2phi,expphi,err]=plotd3d(u,f,solver,argcaxis,mnp,func)

%
% $Id: plotd3d.m,v 1.2 2011/03/26 12:56:38 patrick Exp $
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

x=linspace(solver.xa,solver.xb,solver.nx);
y=linspace(solver.yc,solver.yd,solver.ny);
z=linspace(solver.ze,solver.zf,solver.nz);

hx=(solver.xb-solver.xa)/(solver.nx-1);
hy=(solver.yd-solver.yc)/(solver.ny-1);
hz=(solver.zf-solver.ze)/(solver.nz-1);

m=mnp(1); n=mnp(2); p=mnp(3);

U=zeros(solver.nx+2,solver.ny+2,solver.nz+2);
U(2:end-1,2:end-1,2:end-1)=u;
U([1 end],2:end-1,2:end-1)=u([1 end],:,:);
U(2:end-1,[1 end],2:end-1)=u(:,[1 end],:);
U(2:end-1,2:end-1,[1 end])=u(:,:,[1 end]);
U([1 end],[1 end],[1 end])=u([1 end],[1 end],[1 end]);
L=6.0*del2(U,hy,hx,hz);

d2phi=L(2:end-1,2:end-1,2:end-1);
d2phi=d2phi(2:end-1,2:end-1,2:end-1);

expphi=-nlF(u,Ne3d(solver,solver.ngrid),Te3d(solver,solver.ngrid));
expphi=expphi(2:end-1,2:end-1,2:end-1);

err=L(2:end-1,2:end-1,2:end-1)+ ...
	nlF(u,Ne3d(solver,solver.ngrid),Te3d(solver,solver.ngrid))-f;
err=abs(err(2:end-1,2:end-1,2:end-1));

x=x(2:end-1);
y=y(2:end-1);
z=z(2:end-1);

if 1

subplot(m,n,p), 
errvis=err;
errvis(find(errvis<eps))=eps;
imagesc(x,y,log10(errvis(:,:,fix(end/2)))'), axis xy
xlabel('x')
ylabel('y')
colorbar('v');
if ~isempty(argcaxis),
  caxis(argcaxis); colorbar
end
if mod(p-1,n)~=0, set(gca,'yticklabel',[]); end;
%if fix((p-1)/n)~=m-1, set(gca,'xticklabel',[]); end;
title(sprintf('log10 |d(\\Phi%s)|',func));

else

subplot(m,n,p), imagesc(x,y,f(:,:,fix(end/2))), colorbar('h');
if mod(p-1,n)~=0, set(gca,'yticklabel',[]); end;
if mod(p-1,n)==0, ylabel(func); end;
if fix((p-1)/n)~=m-1, set(gca,'xticklabel',[]); end;
if fix((p-1)/n)==0, title('\rho'); end;

subplot(m,n,p+1), imagesc(x,y,d2phi(:,:,fix(end/2))), colorbar('h');
if mod(p,n)~=0, set(gca,'yticklabel',[]); end;
if fix(p/n)~=m-1, set(gca,'xticklabel',[]); end;
if fix(p/n)==0, title('\nabla^2\Phi'); end;

subplot(m,n,p+2), imagesc(x,y,-expphi(:,:,fix(end/2))), colorbar('h');
if mod(p+1,n)~=0, set(gca,'yticklabel',[]); end;
if fix((p+1)/n)~=m-1, set(gca,'xticklabel',[]); end;
if fix((p+1)/n)==0, title('-exp(\Phi/T)'); end;

subplot(m,n,p+3), imagesc(x,y,err(:,:,fix(end/2))), colorbar('h');
if mod(p+2,n)~=0, set(gca,'yticklabel',[]); end;
if fix((p+2)/n)~=m-1, set(gca,'xticklabel',[]); end;
if fix((p+2)/n)==0, title('\nabla^2\Phi-exp(\Phi/T)+\rho'); end

end

drawnow

fprintf(1,'defect(%s) m %.2e std %.2e\n',func,mean(err(:)),std(err(:)));

