% Calculate image sample area assuming Heymsfield and Parish (1978)
%   bins_mid - mid-point of each bins in doide number
%   res - photodiode resolution, bin width in microns 
%   armdst - distance between probe arms in millimeters
%   num_diodes - number of photodiodes (does not need to equal number of bins)
%   SAmethod - method to calculate SA. Could be 0: center in, 1: entire in and 2 with correction
%       
% ** Created to replace calc_sa to include different bin setup, and three choices of SA
%    calculation methods. Notice the parameter differences.        Will, 2014/06/04
function sa = calc_sa_randombins(bins_mid,res,armdst,num_diodes, SAmethod, probetype)

% calculates OAP SA in mm^2 with res in um, armdist in mm


res = res * 1e-3; %mm
%max_diameter = max_diameter*.5e-3;


radius = bins_mid/2; %radius = bins_mid .* res/2;  You can use this one if you provide midpoint in doide number 
diameter = 2 * radius;

% Calculate the width
switch SAmethod
    case 0
        % Center in 
        EAWci = num_diodes*res;
        EAWri = EAWci;
    case 1
        % Entire in
        EAWci = (num_diodes-(bins_mid/res)+1)*res;
        EAWci(EAWci<0)=0;
        EAWri = EAWci;
    case 2
        % With Correction
        EAWci = num_diodes*res;
        EAWri = EAWci + 0.72*diameter; 
end

% Calculate the DOF
lambda = 680 * 1e-6;             % mm,laser wavelength
if probetype==2
    DOF =20.52*radius.^2*1000; 
else
    DOF = 6*radius.^2/lambda;        % Using 
end
    
DOF(DOF > armdst) = armdst;

sa = DOF .* EAWri;
