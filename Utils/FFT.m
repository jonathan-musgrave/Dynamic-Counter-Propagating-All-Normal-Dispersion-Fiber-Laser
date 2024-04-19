% FFT function for a step that will keep all units correct and so you don't
% have to think about it
function uu = FFT(UU, nt, dt, fRep)
    % Note that it is shifted back because the frequency domain is
    % specified as a shifted axis
    % This is done in order to ensure that the sum(abs(UU).^2).*dt == PowerAvg
    % at each dz step
    uu = fftshift(fft((UU)).*sqrt(1./(nt*dt*fRep)));
    

end