function[Frot]=rotavg0(Fq,th,ph)

%========================================
% Calc. rot-avg,
Iph=sum(Fq,3);                     % sum over phi
for i=1:length(th)
    Iph(:,i)=Iph(:,i)*sin(th(i));  
end

dth=th(2)-th(1);
dph=2*pi/length(ph); % for some reason this is slightly different to ph(2)-ph(1)

% sum over theta, multiply by dth*dph, divide by 4*pi
Frot=sum(Iph,2)*dth*dph*.25*pi^-1; 
Frot=Frot';
%========================================

return