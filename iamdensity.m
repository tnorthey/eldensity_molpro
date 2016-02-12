function [] = iamdensity2(mldfile,dx,xMin,xMax,yMin,yMax,zMin,zMax)
%IAMDENSITY Reads .xyz file of molecule, loads spline functions
% of each atomic radial electron density. 
% Creates IAM molecular density on a 3D grid of size N^3.

addpath('myfunctions')

[~,title0,~]=fileparts(mldfile);
title_=strcat(title0,'_elDensity_iam');

[~,~,~,~,~,~,~,~,~,~,~,~,~,Atoms] = mldread2(mldfile);
% Atoms=xyzread(xyzfile);
disp(Atoms)
% au2ang = 0.52917721092d0;
coords=Atoms(:,3:5);
ZZ=Atoms(:,2);
Nat=size(Atoms,1);

% keyboard

load spline_atomic_radial

% Create a grid
xGrid=xMin:dx:xMax;
yGrid=yMin:dx:yMax;
zGrid=zMin:dx:zMax;

[X,Y,Z] = meshgrid(xGrid,yGrid,zGrid);

Nx=length(xGrid);
Ny=length(yGrid);
Nz=length(zGrid);

Rho3D=zeros(Nx,Ny,Nz);
Rho3D=squeeze(Rho3D);
for a=1:Nat
    disp(100*a/Nat)
    r=sqrt((X-coords(a,1)).^2+(Y-coords(a,2)).^2+(Z-coords(a,3)).^2);
    r(r>6)=6.0;  % > 6 is outside interpolant range
    if ZZ(a)==1
        Rho_=ppval(pH,r);
    elseif ZZ(a)==6
        Rho_=ppval(pC,r);
    elseif ZZ(a)==8
        Rho_=ppval(pO,r);        
    end
    Rho3D=Rho3D+squeeze(Rho_);
end
clear Rho_

save(['results/',title_],'Rho3D','xGrid','yGrid','zGrid','Atoms')

% Total numer of electrons; dx^3 is the volume element
totRho=sum(sum(sum(Rho3D)))*dx^3;
disp('Integrated number of electrons: ');
disp(num2str(totRho));
disp('True number of electrons: ');
disp(num2str(sum(Atoms(:,2))));
% 
% keyboard
% 
% % plot
% contourf(zGrid,yGrid,log10(squeeze(Rho3D(:,ceil(Nz/2),:))),100,'edgecolor','none');
% colorbar;
% xlabel('$y (a_0)$','interpreter','latex');
% ylabel('$x (a_0)$','interpreter','latex');
% axis equal
% axis tight
% 
% keyboard

end
