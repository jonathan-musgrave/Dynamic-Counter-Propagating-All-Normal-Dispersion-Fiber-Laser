function UU = IFFT(uu, nt, dt,fRep)
    % Note that it is not shifted back because the frequency domain is
    % specified as a shifted axis
    % This is done in order to ensure that the sum(abs(UU).^2) == PowerAvg
    % at each point
    UU = ifft(ifftshift(uu)).*sqrt(dt.*fRep.*nt);
end