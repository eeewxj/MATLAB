function [Eout,gain] = amp_simp(Ein,GssdB,PoutsatdB,NF)
% amp_simp(Ein,GssdB,PoutsatdB,NF)
% simple model of optical amplifier. The model includes the gain saturation and the ASE noise

% amplifier parameters:
%   small signal gain: GssdB (dB)
%   output saturation power:PoutsatdB (dBm)
%   noise figure: NF (dB)

% Simulation parameters:
%   Gain tolerance: tol, used as the threshold value to exit the gain calculation loop

% The input is a column vector containing block N samples of the optical
% signal sampling at the rate 1/Ts
% The output is calculated using
%   Eout = Ein*sqrt(G) + Enoise
% where: G is the saturation gain
%        G = Gss*exp(-(G-1)*Pin/Psat) (eq1)
% Enoise id the complex noise with the noise power
%       Pase = (10^(NF/10))*(G-1)hf/2*BW;
%       BW = 1/Ts;

global Ts;      % sampling period
global Fcar;
h_plan = 6.626e-34;

Gss = 10^(GssdB/10);
Poutsat = (10^(PoutsatdB/10))/1000;
Psat = Poutsat*(Gss-2)/Gss/log(2);
% Pinsat = 2*Poutsat/Gss;

N = size(Ein,1);
Pin = mean(Ein.*conj(Ein));

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
gain = G;
Egain = sqrt(G)*Ein;

BW = 1/Ts;
FigNoise = 10^(NF/10);
Pase = (FigNoise*G-1)*h_plan*Fcar/2*BW;
PasedB = 10*log10(Pase);

Eout = Egain + wgn(N,1,PasedB,'complex');
