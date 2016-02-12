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

for k=2:nr                     % loop through spheres of non-zero radius r(kk)
    for j=1:nth                % theta loop 
        for i=1:nph            % phi loop, note: skips 2*pi as f(0)=f(2*pi)       
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