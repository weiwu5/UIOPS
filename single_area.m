function area=single_area(diameter, habit)
% 
% Used to calculate the particle area using A-D relations
% This is in cgs units
%  
% Created on Feb 14, 2014 by Will Wu
%
area = 0;
switch char(habit)
    %case 't'  % Sphere and Tiny
        %iwc=0.91*3.1415926/6*diameter^3;
        %a=0.049;
        %b=2.8;
        %iwc = a *diameter^b;
    case {'l','o'}  % Linear and Oriented
        if diameter<0.03
            a=0.0696;
            b=1.5;
        else
            a=0.0512; %0.000907;
            b=1.414;
        end
        area = a *diameter^b;
    case 'h'  % Plate
        a=0.65;
        b=2.0;
        area = a *diameter^b;
    case {'i','a','s','t'}       % Inregular 
        a=0.2285;
        b=1.88;
        area = a *diameter^b;
    case {'g'}       % Graupel 
        a=0.5;
        b=2.0;
        area = a *diameter^b;
    case 'd'       % Dendrite  
        a=0.21;
        b=1.76;
        area = a *diameter^b;
%     case 'a'       % Aggregate 
%         a=0.0033;
%         b=2.2;
%         iwc = a *diameter^b;
end

area = area*100; % Change unit into mm^2

end