function [NW,NSpec] = OSA_Resol(w1,w2,OSpec,resol,lamda0)

% wavelength unit is nm, and lamda0 unit is m

global v c_const;

NW = (w1:((w2-w1)/1000):w2);  % nm
a = size(OSpec,1);

NSpec = zeros(a,1001);
for ii = 1:1001
    f3dB = resol*c_const*1e9/NW(ii)^2;       % Hz
    
    deltavN = c_const/lamda0-c_const/NW(ii)/1e-9;        % Hz
    [~,aux1b] = min(abs(v - deltavN));
    
    Spec = circshift(OSpec',-aux1b);
    for pp = 1:a
        Spec1 = filter_gaus(fftshift(Spec(:,pp)),f3dB,2);
        NSpec(pp,ii) = 10*log10(sum(Spec1));
    end
end