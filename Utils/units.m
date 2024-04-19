

% km = 1e12; m = 1e9; cm = 1e7; mm = 1e6; um = 1e3; nm = 1; pm =1e-3; 
% fs = 1e-3; ps = 1; ns = 1e3; us = 1e6; ms = 1e9; s = 1e12; 
% Hz = 1./s; KHz = 1./ms; MHz = 1/us; GHz = 1/ns; THz = 1/ps; 
% nW = 1e-9; uW = 1e-6; mW = 1e-3; W = 1; kW = 1e3; MW = 1e6; GW = 1e9;
% nJ = nW*s; uJ = uW*s; mJ = mW*s; J = W*s; KJ = kW*s; MJ = MW*s; GJ = GW*s;



 m = 1; km = 1e3; cm = 1e-2; mm = 1e-3; um = 1e-6; nm = 1e-9; pm =1e-12; 
Hz = 1; KHz = 1e3; MHz = 1e6; GHz = 1e9; THz = 1e12; 
nW = 1e-9; uW = 1e-6; mW = 1e-3; W = 1; 
kW = 1e3; MW = 1e6; GW = 1e9; PW = 1e12; TW = 1e15;
pJ = 1e-12; nJ = 1e-9; uJ = 1e-6; mJ = 1e-3; J = 1; KJ = 1e3; MJ = 1e6; GJ = 1e9;
fs = 1e-15; ps = 1e-12; ns = 1e-9; us = 1e-6; ms = 1e-3; s = 1;  F = 1; kg = 1; K= 1;
V = 1;
%%%%%%%%%%%%% Fundamental Constants %%%%%%%%%%%%%%%%%%%%
eps0 = 8.85418e-12*s^4./(m^3); %s^4A^3/(kg m^3)
mu0 = 1.25663706e-6./m;  %H/m
c = 1./sqrt(eps0*mu0);
h = 6.62607015e-34*m^2; hbar = h/(2*pi); 
deg = pi/180; rad = 1; 
dBm = @(x) 10^(x/10)/1e3; Perc2dB = @(x) 10*log10(x); W2dBm = @(x) 10*log10(x./mW);