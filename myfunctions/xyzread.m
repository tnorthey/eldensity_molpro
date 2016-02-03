function [Atoms] = xyzread(xyzfile)

au2ang = 0.52917721092d0; % convert Ang coordinates to au at the end

disp(['xyzread: Read from file ' xyzfile])
fid = fopen(xyzfile, 'r');                 % Open file for reading
tline = fgetl(fid);                        % Read 1st line from file, removing newline characters
Nat=str2double(tline);

for i=1:2; tline = fgetl(fid); end 

% Atoms:
Atoms=zeros(Nat,5);
for i=1:Nat
    p = textscan(tline,'%s %f %f %f');  
    Atoms(i,:) = [i ele2Z(p{1}) ... 
        p{2}/au2ang...  % Atom Number, no of electrons, x, y, z (x,y,z in Bohr)
        p{3}/au2ang...
        p{4}/au2ang]; 
    tline = fgetl(fid);                    % next line
end 

end % END FUNCTION xyzread

function[Z]=ele2Z(ele)
if strcmpi(ele,'H');
    Z=1;
elseif strcmpi(ele,'He');
    Z=2;
elseif strcmpi(ele,'Li');
    Z=3;
elseif strcmpi(ele,'Be');   
    Z=4;
elseif strcmpi(ele,'B');   
    Z=5;
elseif strcmpi(ele,'C');
    Z=6;
elseif strcmpi(ele,'N');
    Z=7;
elseif strcmpi(ele,'O');
    Z=8;  
elseif strcmpi(ele,'F');
    Z=9; 
elseif strcmpi(ele,'Ne');
    Z=10; 
elseif strcmpi(ele,'S');
    Z=16;    
else
    disp('Element not listed, easy to add.')
end
end % end function
