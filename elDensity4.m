function [ Rho3D, xGrid, yGrid, zGrid ] = elDensity4(mldFileName,morange,dx,xMin,xMax,yMin,yMax,zMin,zMax)

%elDensity3 Calcualte electron density from a file mldFileName
% The function adopts the notation used in the function mldread to calcuate
% the electron denisty for a grid specified by the input arguments
% ============================================================
% Author: Nikola Zotev
% Last Modified: 16/06/2014
% Thomas Northey edit 03/02/2016:
% Added automatic boxsize, variable amount of input arguments,
% auto title, and save to results folder
% ============================================================
%   Input:
% dx    - double - the length of the grid
% xMin  - double - the min value of the grid in x
% xMax  - double - the max value of the grid in x
% yMin  - double - the min value of the grid in y
% yMax  - double - the max value of the grid in y
% zMin  - double - the min value of the grid in z
% zMax  - double - the max value of the grid in z
% ========================================================
%   Output:
% Rho3D - 3D matrix     - the electron density on a Cartesian grid
% xGrid - 1D col vector - the x values of the grid
% yGrid - 1D col vector - the y values of the grid
% zGrid - 1D col vector - the z values of the grid
% =======================================================
% Future modifications:
% 1) the function is based on an old version of mldread() => update
% 2) !!! plot doesn't work if xMin!=yMin and xMax!=yMax 
% =====================================================


% Read from mldread
% Col vectors:
% bi    = occupancy
% M     = Mol. orb. expansion coefficient
% ga    = exponents of primitive Gaussian
% ci    = coefficinet in front of primitve Gaussian
% l,m,n = components of angular momentum
% xx, yy, zz = position of atoms == centre of Gaussians
% ppmo  = primitives per MO
% moocc = number of MO
% primlist = not used ???
% Atoms = matrix for all atoms in format [Atom Number, no of electrons, x, y, z]

% EDIT by TN, add title for saving results later,
% and add speed up if statement, omit calculation of grid points far away
% (>20 au) from centre of primitive.

% switch nargin
%     case 3
%         mldFileName=varargin{1};
%         Nx=2*ceil(varargin{3}/2.0)-1; % make odd so the centre of an atom at the origin is exactly the central point on the grid
%         method=varargin{2};
%     case 2
%         mldFileName=varargin{1};
%         Nx=79; 
%         method=varargin{2};
%     case 1
%         mldFileName=varargin{1};
%         Nx=79;  % good value, tested for H and C atoms
%         method='ai';
%     otherwise
%         print Need at least 1 input: mldfile 
% end
% Ny=Nx; Nz=Nx; %xMin, xMax, yiMin, yMax, zMin, zMax

addpath('myfunctions')

[~,title0,~]=fileparts(mldFileName);
title_=strcat(title0,'_elDensity');

[ bi,M,ga,ci,l,m,n,xx,yy,zz,ppmo,moocc,~,Atoms] = mldread2(mldFileName);

% xyz=Atoms(:,3:5);
% ZZ=Atoms(:,2);
% Nat=size(Atoms,1);

% show Atoms
disp(Atoms)

% Make automatic cubic box limits
% Haven't really tested to failure (if any):
% lxyz=size(xyz,1);
% dist=zeros(lxyz^2,1); dist0=dist;
% cc=0;
% for i=1:lxyz
%     for j=1:lxyz
%         cc=cc+1;
%         dist0(cc) = sqrt(sum(xyz(i,:).^2));  % square the row, sum, then sqrt
%         dist(cc)  = sqrt(sum((xyz(i,:)-xyz(j,:)).^2));
%     end
% end
% max0=max(dist0);   % furthest distance from the origin
% maxdist=max(dist); % longest atom-atom distance
% atom_radius=7.0;   % gives integral ~99.5% to true Nelec for H and C atoms
% 
% if max0>maxdist                  % then the origin is outside the molecule
%     boxsize=max0+atom_radius;    % and my symmetric about the origin boxsize is defined by distance from the origin
% elseif maxdist>max0              % the origin is inside the molecule
%     boxsize=maxdist+atom_radius;
% else                             % there's only 1 atom at the origin
%     boxsize=atom_radius;         
% end
% disp(strcat('Cubic boxsize (bohr), symmetric about origin: ',num2str(boxsize)))
% xMin=-0.5*boxsize; yMin=xMin; zMin=xMin;  
% xMax= 0.5*boxsize; yMax=xMax; zMax=xMax; 

% Create a grid
xGrid=xMin:dx:xMax;
yGrid=yMin:dx:yMax;
zGrid=zMin:dx:zMax;

[X,Y,Z]=meshgrid(xGrid,yGrid,zGrid);

% Reserve memory
Nx=length(xGrid);
Ny=length(yGrid);
Nz=length(zGrid);
Rho3D=zeros(Nx,Ny,Nz);

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
contourf(zGrid,yGrid,squeeze(Rho3D(ceil(Nz/2),:,:)),50,'edgecolor','none');
colorbar;
xlabel('$y (a_0)$','interpreter','latex');
ylabel('$x (a_0)$','interpreter','latex');
axis equal
axis tight

keyboard

end



