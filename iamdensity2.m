function [] = iamdensity2(mldfile,dx,xMin,xMax,yMin,yMax,zMin,zMax)
%IAMDENSITY Reads .xyz file of molecule, loads spline functions
% of each atomic radial electron density. 
% Creates IAM molecular density on a 3D grid of size N^3.

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

load c_6-31g_hf_elDensity
RhoC=Rho3D;
xC=xGrid; yC=yGrid; zC=zGrid;
load h_6-31g_hf_elDensity
RhoH=Rho3D;
xH=xGrid; yH=yGrid; zH=zGrid;
load o_6-31g_hf_elDensity
RhoO=Rho3D;
xO=xGrid; yO=yGrid; zO=zGrid;

% Create a grid
xGrid=xMin:dx:xMax;
yGrid=yMin:dx:yMax;
zGrid=zMin:dx:zMax;

[X,Y,Z] = meshgrid(xGrid,yGrid,zGrid);

Nx=length(xGrid);
Ny=length(yGrid);
Nz=length(zGrid);

Rho3D=zeros(Nx,Ny,Nz);
for a=1:Nat
    disp(100*a/Nat)
    X_=X+coords(a,1); Y_=Y+coords(a,2); Z_=Z+coords(a,3); 
    if ZZ(a)==1
        Rho_=interp3(xH,yH,zH,RhoH,X_,Y_,Z_);
        Rho_(isnan(Rho_))=0;
        Rho3D=Rho3D+Rho_;    
    elseif ZZ(a)==6
        Rho_=interp3(xC,yC,zC,RhoC,X_,Y_,Z_);
        Rho_(isnan(Rho_))=0;
        Rho3D=Rho3D+Rho_; 
    elseif ZZ(a)==8
        Rho_=interp3(xO,yO,zO,RhoO,X_,Y_,Z_); 
        Rho_(isnan(Rho_))=0;
        Rho3D=Rho3D+Rho_; 
    end
end

save(['results/',title_],'Rho3D','xGrid','yGrid','zGrid','Atoms')

% Total numer of electrons; dx^3 is the volume element
dx=abs(xGrid(2)-xGrid(1));
totRho=sum(sum(sum(Rho3D)))*dx^3;
disp('Integrated number of electrons: ');
disp(num2str(totRho));
disp('True number of electrons: ');
disp(num2str(sum(Atoms(:,2))));

% plot
contourf(zGrid,yGrid,squeeze(Rho3D(ceil(Nz/2),:,:)),100,'edgecolor','none');
colorbar;
xlabel('$y (a_0)$','interpreter','latex');
ylabel('$x (a_0)$','interpreter','latex');
axis equal
axis tight

keyboard

end
