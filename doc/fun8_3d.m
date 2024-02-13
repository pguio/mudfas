function [u,f,solver]=fun8(solver)
% function [u,f,solver]=fun8(solver)

%
% $Id: fun8_3d.m,v 1.2 2011/03/26 12:56:38 patrick Exp $
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

global iprer ipost

solver.iprer=iprer;
solver.ipost=ipost;

u=zeros(solver.nx,solver.ny,solver.nz);

global a kx ky kz X Y Z phit

x=linspace(0,1,solver.nx);
y=linspace(0,1,solver.ny);
z=linspace(0,1,solver.nz);
[X,Y,Z]=ndgrid(x,y,z);

nb=1.2;
f=ones(solver.nx,solver.ny,solver.nz);
imn=fix(solver.nx/3);
imx=solver.nx-imn;
jmn=fix(solver.ny/2);
jmx=solver.ny;
kmn=fix(solver.nz/3);
kmx=solver.nz-kmn;

if 0
f(imn:imx,jmn:jmx,kmn:kmx)=nb*f(imn:imx,jmn:jmx,kmn:kmx);
else
xc=0.5*(x(imn)+x(imx));
zc=0.5*(z(kmn)+z(kmx));
D=sqrt((X-xc).^2+(Z-zc).^2);
II=find(D<=.5*(x(imx)-x(imn)) & Y>=y(jmn));
f(II)=nb*f(II);
end

randn('state',0);
m=0.0;
sd=0.0;
f=f+sd*(randn(solver.nx,solver.ny,solver.nz)-m);

phit=zeros(size(f));
phit(imn:imx,jmn:jmx,kmn:kmx)=solver.temin*log(nb);

