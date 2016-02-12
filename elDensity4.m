function [ Rho3D, xGrid, yGrid, zGrid ] = elDensity4(mldFileName,morange,dx,xMin,xMax,yMin,yMax,zMin,zMax,sph_on)

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

% Cartesian or Spherical 
if sph_on==0
    [X,Y,Z]=meshgrid(xGrid,yGrid,zGrid);
elseif sph_on==1
    rmax=max(xGrid);   % not so general...
    C=sphcoords(Nx,Ny,Nz,rmax);
    X=C{1}; Y=C{2}; Z=C{3};
end

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
if sph_on==0
    dx=abs(xGrid(2)-xGrid(1));
    totRho=sum(sum(sum(Rho3D)))*dx^3;
elseif sph_on==1
    r=linspace(0,rmax,Nx);
    dr=r(2)-r(1);
    th_=linspace(0,pi,Ny);
    dth=th_(2)-th_(1);
    th=[];
    for i=1:Nx
        th=[th;th_];
    end
    ph=linspace(0,2*pi,Nz);
    dph=ph(2)-ph(1);
    
    totRho=sum(sum(sum(Rho3D,3).*sin(th),2).*r'.^2,1)*dr*dth*dph;
end
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


function[C]=sphcoords(Nx,Ny,Nz,rmax)
%=======================================================================
% spherical coordinate system
nr=Nx;
nth=Ny;
nph=Nz;

r=linspace(0,rmax,nr);       % radius
th=linspace(0,pi,nth);       % theta
ph=linspace(0,2*pi,nph);     % phi

qx=zeros(nr,nth,nph);qy=qx;qz=qx;

for k=1:nr                     % loop through spheres of non-zero radius r(kk)    
    for i=1:nph                % phi loop, note: skips 2*pi as f(0)=f(2*pi)
        for j=1:nth            % theta loop           
            % Create x,y,z as a function of spherical coords...           
            qx(k,j,i)=r(k)*sin(th(j))*cos(ph(i));
            qy(k,j,i)=r(k)*sin(th(j))*sin(ph(i));
            qz(k,j,i)=r(k)*cos(th(j));        
        end
    end  
end
C=cell(3,1); C{1}=qx; C{2}=qy; C{3}=qz; % coordinate matrix
%=======================================================================
return

