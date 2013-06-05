function gain = gain_saturated(Pin,GssdB,PoutsatdB)
% calculation the gain of the amplifier given the input power, the small signal gain and saturation power

% where: G is the saturation gain
%        G = Gss*exp(-(G-1)Pin/Psat) (eq1)

% Amplifier parametes:
%   GssdB = 20          % (dB)
%   PoutsatdB = 10      % (dBm)
%   NF = 7              % (dB)

Gss = 10^(GssdB/10);
Poutsat = (10^(PoutsatdB/10))/1000;
Psat = Poutsat*(Gss-2)/Gss/log(2);
% Pinsat = 2*Poutsat/Gss;

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