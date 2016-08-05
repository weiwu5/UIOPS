%% Returns terminal velocity for a single particle
%  Both options to calculate the terminal velocity
%  Default is to use the Heymsfield and Westbrook (2010) method, 
%  but you can also choose to to use Mitchel (1996)
%  Created by Will Wu, 2014/01/15
%  - Mass and Diameter uses metric system
%  - Pressure use hPa
%  - Temperature use Celsius
function vt = single_vt(diameter, area_ratio, mass, P, T)

usingMithcell=0; % Setting 0 to use Heymsfield method, other value for Mitchell method
g=9.8;
pi=3.1415926;

% Calculate environmental conditions
T = T + 273.15;
P = P*100;
rho_a = P/(287.15*T);
eta = 18.27*(291.15+120)./(T+120)*(T/291.15)^(3/2)/10^(6);  % Sutherland's formula to calculate the dynamic viscosity
nu = eta/rho_a;  % kinectic viscosity

if 0==usingMithcell 
    % Calculate modified Best Number using Heymsfield and Westbrook (2010).
    % using drag C=0.35 and epsilon=8.0
    X=rho_a/(eta^2)*8*mass*g/(pi*area_ratio^0.5);
    ReynoldN=16*(sqrt(1+4*sqrt(X)/64/sqrt(0.35))-1)^2;
else 
    % Calculate modified Best Number using Mitchell (1996). 
    % This is actually a special case of Heymsfield and Westbrook (2010), with
    % k=0, and here we use drag C=0.6 and epsilon=5.83
    % We calculate from the original equations without using power law
    % approximation
    X=rho_a/(eta^2)*8*mass*g/(pi*area_ratio);
    ReynoldN=5.83^2/4*(sqrt(1+4*sqrt(X)/(5.83^2)/sqrt(0.6))-1)^2;
end

vt=nu/diameter*ReynoldN;

end

