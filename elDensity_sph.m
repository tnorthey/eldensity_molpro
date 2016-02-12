function [Rho3D,r,th,ph] = elDensity4_sph(mldFileName,morange,dr,rMax)

addpath('myfunctions')

[~,title0,~]=fileparts(mldFileName);
title_=strcat(title0,'_elDensity_sph');

[ bi,M,ga,ci,l,m,n,xx,yy,zz,ppmo,moocc,~,Atoms] = mldread2(mldFileName);

% show Atoms
disp(Atoms)

% Create a grid
r=0:dr:rMax;

% Reserve memory
Nr=length(r);
Rho3D=zeros(Nr,Nr,Nr);

% polar and azimuthal
th=linspace(0,pi,Nr);
ph=linspace(0,2*pi,Nr);

% spherical coordinates,
C=sphcoords(Nr,Nr,Nr,rMax);
X=C{1}; Y=C{2}; Z=C{3};

% Normalisation factor in front of each prim. Gaussian
normA=(2/pi)^0.75;
normB=2.^(l+m+n);
normC=ga.^((2*l+2*m+2*n+3)/4);
normD=factd(2*l-ones(size(l))) .* factd(2*m-ones(size(l))) .* factd(2*n-ones(size(l)));

norm=normA .* normB .* normC ./ normD.^ 0.5;

% Sum primitives that belong to a given MO
for mo=morange
    
    disp(100*mo/length(morange))
    
    Psi=zeros(Nr,Nr,Nr);
    
    for i=(mo-1)*ppmo+1:mo*ppmo

        r_sq=(X-xx(i)).^2+(Y-yy(i)).^2+(Z-zz(i)).^2;

        Psi = Psi + norm(i)*M(i)*sqrt(bi(i))*ci(i).*...
            (X-xx(i)).^l(i).*(Y-yy(i)).^m(i).*(Z-zz(i)).^n(i).*...
            exp(-ga(i).*r_sq);

    end    
    % Square to get density
    Rho3D=Rho3D+Psi.^2;
end

% check total electrons,
th_=zeros(Nr,Nr);
for i=1:Nr
    th_(i,:)=th;
end
totRho=sum(sum(sum(Rho3D,3).*sin(th_),2).*r'.^2,1)*rMax*2*pi^2/(Nr*(Nr-1)^2);
disp('Integrated number of electrons: ');
disp(num2str(totRho));
disp('True number of electrons: ');
disp(num2str(sum(Atoms(:,2))));

% keyboard

% rotationally average,
Rhor=rotavg0(Rho3D,th,ph);

% save file in results folder,
title0=strcat('results/',title_);
disp('saving .mat file:')
disp(title0)
save(title0,'r','Rhor')
    
return