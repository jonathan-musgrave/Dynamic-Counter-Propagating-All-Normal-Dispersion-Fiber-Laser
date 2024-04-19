% returns C.O.M of an array axis given the dependent array

% Useful for calculating center frequency or pulse time


function [xCOM] = calc_COM(x,y)
    tot_mass = sum(y(:));
    xCOM = sum(x.*y)./tot_mass;
end
