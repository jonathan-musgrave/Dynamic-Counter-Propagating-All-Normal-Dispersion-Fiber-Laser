
% This function interpolates and collectects the YDF Spectroscopy data for
% the fiber definition over the frequency axis that is given as an input.
function [emi_sig, abs_sig] = getYDFSpectroData_v2(W)
    % addpath '../Utils/YDF Spectroscopic Data'
    Data = readmatrix('YDFSpectroData_2.xls');
    YbcsAbsLam = Data(:,1); % in m
    YbcsAbsSig = Data(:,3); % m^2
    YbcsEmiLam = Data(:,1);
    YbcsEmiSig = Data(:,4); 
    
    YbcsEmiLam = [YbcsEmiLam;1.151e-6]; % Adding some extra tails for interpolation
    YbcsEmiSig = [YbcsEmiSig;0];
    
    YbcsAbsLam = [.849e-6;YbcsAbsLam]; % Adding some extra tails for interpolation
    YbcsAbsSig = [0;YbcsAbsSig];
    c = 299792458;
    abs_sig = interp1(2*pi*c./YbcsAbsLam,YbcsAbsSig,W);
    abs_sig(isnan(abs_sig)) = 0;
    emi_sig = interp1(2*pi*c./YbcsEmiLam,YbcsEmiSig,W);
    emi_sig(isnan(emi_sig)) = 0;

end