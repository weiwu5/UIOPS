function iwc=single_mass(diameter, habit)
% This is in cgs units
iwc = 0;
switch char(habit)
    case {'s','t'}  % Sphere and Tiny
        %iwc=0.91*3.1415926/6*diameter^3;
        a=0.00294;
        b=1.9;
        %a=0.049;
        %b=2.8;
        iwc = a *diameter^b;
    case {'l','o'}  % Linear and Oriented
        if diameter<0.03
            a=0.001666;
            b=1.91;
        else
            a=0.000907;
            b=1.74;
        end
        iwc = a *diameter^b;
    case 'h'  % Plate
        a=0.00739;
        b=2.45;
        iwc = a *diameter^b;
    case {'i','a'}       % Inregular 
        a=0.00294;
        b=1.9;
        iwc = a *diameter^b;
    case 'g'       % Graupel 
        a=0.049;
        b=2.8;
        iwc = a *diameter^b;
    case 'd'       % Dendrite  
        a=0.000516;
        b=1.8;
        iwc = a *diameter^b;
%     case 'a'       % Aggregate 
%         a=0.0033;
%         b=2.2;
%         iwc = a *diameter^b;
end

end
