function [slotfreq,powerinres] = resolution(u0,vs,res,nt,dt)

% This function gives the total (sum) of u0 in a badwidth res.
% It also gives the centerfrequency of the interval at which
% the sum was done(slotfreq). The frequency 0 is always forced
% to appear in alotfreq and the other frequencies are at steps
% +res and -res from 0. vs are the frequencies of fft(u0).
% 
% INPUTS
% 
% u0 - initial field (W^.5)
% vs - frequencies of fft(u0) in Ghz
% res - the required resolution bandwidth (Ghz)
% nt - length(u0)
% dt - distance between two  sucesive points
% 
% OUTPUTS
% 
% slotfreq - center frequencies of each interval in which the sum is done (GHz)
% powerinres - total power (fft(u0))^2 of each interval with bandwidth res (W)

pointsperslot = res*nt*dt/1000;
spectralpower = fftshift((abs(fft(u0)).^2)/nt^2);
slotfreq = [(0:res:max(vs)),(-res:-res:min(vs))]; % vector of center freq of each slot (GHz)
slotfreq = sort(slotfreq);

if min(vs) < min(slotfreq)-res/2
startcountat = max(find(vs < (min(slotfreq)-res/2)));
else
    aux1=[];    
    for i=2:length(slotfreq)-1,    
    aux1(i-1) = slotfreq(i);    
    end
    slotfreq = aux1;
    startcountat = max(find(vs < (min(slotfreq)-res/2)));
end

powerinres=[];
aux=0;
for j = 0:length(slotfreq)-1,
   for i =1:floor(pointsperslot),
   aux = aux + spectralpower(startcountat + floor(j*pointsperslot) + i);
   end
 aux = aux + (pointsperslot - floor(pointsperslot))*spectralpower(startcountat + floor(j*pointsperslot));
 powerinres(j+1) = aux;
 aux=0;
end