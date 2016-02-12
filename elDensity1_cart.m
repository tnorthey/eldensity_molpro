function [ Rho3D, xGrid, yGrid, zGrid ] = elDensity4_cart(mldFileName,morange,dx,xMin,xMax,yMin,yMax,zMin,zMax)

addpath('myfunctions')

[~,title0,~]=fileparts(mldFileName);
title_=strcat(title0,'_elDensity');

[ bi,M,ga,ci,l,m,n,xx,yy,zz,ppmo,moocc,~,Atoms] = mldread2(mldFileName);

% show Atoms
disp(Atoms)

% Create a grid
xGrid=xMin:dx:xMax;
yGrid=yMin:dx:yMax;
zGrid=zMin:dx:zMax;

% Reserve memory
Nx=length(xGrid);
Ny=length(yGrid);
Nz=length(zGrid);
Rho3D=zeros(Nx,Ny,Nz);

[X,Y,Z]=meshgrid(xGrid,yGrid,zGrid);

% Normalisation factor in front of each prim. Gaussian
normA=(2/pi)^0.75;
normB=2.^(l+m+n);
normC=ga.^((2*l+2*m+2*n+3)/4);
normD=factd(2*l-ones(size(l))) .* factd(2*m-ones(size(l))) .* factd(2*n-ones(size(l)));

norm=normA .* normB .* normC ./ normD.^ 0.5;

% Sum primitives that belong to a given MO
for mo=morange
    
    disp(100*mo/length(morange))
    
    Psi=zeros(Nx,Ny,Nz);
    
    for i=(mo-1)*ppmo+1:mo*ppmo

        r_sq=(X-xx(i)).^2+(Y-yy(i)).^2+(Z-zz(i)).^2;

        Psi = Psi + norm(i)*M(i)*sqrt(bi(i))*ci(i).*...
            (X-xx(i)).^l(i).*(Y-yy(i)).^m(i).*(Z-zz(i)).^n(i).*...
            exp(-ga(i).*r_sq);

    end    
    % Square to get density
    Rho3D=Rho3D+Psi.^2;
end

% save file in results folder,
title0=strcat('results/',title_);
disp('saving .mat file:')
disp(title0)
save(title0,'xGrid','yGrid','zGrid','Rho3D')
    
% Total numer of electrons; dx^3 is the volume element
dx=abs(xGrid(2)-xGrid(1));
totRho=sum(sum(sum(Rho3D)))*dx^3;
disp('Integrated number of electrons: ');
disp(num2str(totRho));
disp('True number of electrons: ');
disp(num2str(sum(Atoms(:,2))));

% plot
% contourf(zGrid,yGrid,squeeze(Rho3D(ceil(Nz/2),:,:)),50,'edgecolor','none');
% colorbar;
% xlabel('$y (a_0)$','interpreter','latex');
% ylabel('$x (a_0)$','interpreter','latex');
% axis equal
% axis tight

% keyboard

return