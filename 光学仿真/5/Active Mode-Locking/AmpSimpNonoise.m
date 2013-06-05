function [Eout,gain] = AmpSimpNonoise(Ein,GssdB,PoutsatdB)
% amp_Simp(Ein,GssdB,PoutsatdB,NF)
% simple model of optical amplifier. The model includes the gain saturation without noise

% Amplifier parameters:
%   small signal gain: GssdB (dB)
%   output saturation power PoutsatdB (dBm)

% The input is a column vector containing block N samples of the optical signal sampling at the rate 1/Ts
% The output is calculated using
%   Eout = Ein*sqrt(G)
% where: G is the saturated gain
%        G = Gss*exp(-(G-1)Pin/Psat) (eq1)

Gss = 10^(GssdB/10);
Poutsat = (10^(PoutsatdB/10))/1000;
Psat = Poutsat*(Gss-2)/Gss/log(2);
% Pinsat = 2*Poutsat/Gss;

N = size(Ein,1);
Pin = (sum(Ein.*conj(Ein))/N);

% numerical calculation of G from the equation G = (Gss - InG)*Psat/Pin + 1
tol = 0.05;     % tolerence for G calculation
step = Gss/2;
G = Gss;
err = 10;
while (err > tol)
    G1 = Gss*exp(-(G-1)*Pin/Psat);
    err = G1 - G;
    if err>0
        if step<0
            step = -step/2;
        end
    else
        if step>0
            step = -step/2;
        end
        err = -err;
    end
    G = G + step;
end
G = G - step;

Eout = sqrt(G)*Ein;
gain = G;