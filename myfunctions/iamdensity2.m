function [] = iamdensity2(title_,xyzfile,Nx,Ny,Nz,xmin,xmax,ymin,ymax,zmin,zmax)
%IAMDENSITY Reads .xyz file of molecule, loads spline functions
% of each atomic radial electron density. 
% Creates IAM molecular density on a 3D grid of size N^3.

Atoms=xyzread(xyzfile);
coords=Atoms(:,3:5);
ZZ=Atoms(:,2);
Nat=size(Atoms,1);

load spline_atomic_radial

xlen=linspace(xmin,xmax,Nx);
ylen=linspace(ymin,ymax,Ny);
zlen=linspace(zmin,zmax,Nz);
vijk=zeros(Nx,Ny,Nz); r=vijk;
for a=1:Nat
    disp(100*a/Nat)
    for i=1:Nx
        for j=1:Ny
            for k=1:Nz
                % probably can avoid this double loop too...
                r(i,j,k)=sqrt(sum(([xlen(i) ylen(j) zlen(k)]-coords(a,:)).^2)); % radius            
            end
        end
    end
    
    if ZZ(a)==1
        vijk=vijk+ppval(r,pH);
    elseif ZZ(a)==6
        vijk=vijk+ppval(r,pC);
    elseif ZZ(a)==7
        vijk=vijk+ppval(r,pN);
    elseif ZZ(a)==8
        vijk=vijk+ppval(r,pO); 
    elseif ZZ(a)==16
        vijk=vijk+ppval(r,pS);                    
    end   
end

save(['results/',title_,'_iamdensity'],'vijk','xlen','ylen','zlen','Atoms')

end