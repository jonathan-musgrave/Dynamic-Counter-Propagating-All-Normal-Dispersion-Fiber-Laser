function [sys, signal, pumps, fiber] = Rate_Equation_NoASEdz(sys, signal, pumps, fiber, dz)
    run('units.m')
    PSDdz1 = [signal.PSD, pumps.Power];
    w = [(sys.w+2*pi*c./(sys.lambda0)),2*pi*c./pumps.lambda];
    
    if dz<0
        N2 = fiber.N2.*fiber.NT;
        k1       = -fiber.Gamma.*(fiber.emi_sig.*N2 - fiber.abs_sig.*(fiber.NT-N2) - fiber.Loss).*...
             PSDdz1;
        k2       = -fiber.Gamma.*(fiber.emi_sig.*N2-fiber.abs_sig.*(fiber.NT-N2) - fiber.Loss).*...
            (PSDdz1+k1.*dz/2);
        k3       = -fiber.Gamma.*(fiber.emi_sig.*N2-fiber.abs_sig.*(fiber.NT-N2) - fiber.Loss).*...
            (PSDdz1+k2.*dz/2);
        k4       = -fiber.Gamma.*(fiber.emi_sig.*N2-fiber.abs_sig.*(fiber.NT-N2) - fiber.Loss).*...
            (PSDdz1+k3.*dz);
        
        PSDdz2 = PSDdz1 + (dz/6*(k1 + 2*k2 + 2*k3 + k4)).*[ones(1,sys.nt), pumps.BKWchannels];
        pumps.Power(pumps.BKWchannels) = PSDdz2(logical([zeros(1,sys.nt),pumps.BKWchannels]));
       
    elseif dz>0
        R12 = sum(fiber.Gamma.*fiber.abs_sig.*PSDdz1./(hbar.*(w).*fiber.aCore^2.*pi));
        R21 = sum(fiber.Gamma.*fiber.emi_sig.*PSDdz1./(hbar.*(w).*fiber.aCore^2.*pi));
        fiber.N2 = R12./(R12+R21+1./fiber.tau); % in percentage
        N2 = fiber.N2.*fiber.NT; % in Atoms
        % if dz<0  update only backward direction pump/signal and if dz>0
        % update only the forward pump/signal

        k1 = fiber.Gamma.*(fiber.emi_sig.*N2 - fiber.abs_sig.*(fiber.NT-N2) - fiber.Loss).*...
        PSDdz1;
        k2       = fiber.Gamma.*(fiber.emi_sig.*N2-fiber.abs_sig.*(fiber.NT-N2) - fiber.Loss).*...
            (PSDdz1+k1.*dz/2);
        k3       = fiber.Gamma.*(fiber.emi_sig.*N2-fiber.abs_sig.*(fiber.NT-N2) - fiber.Loss).*...
            (PSDdz1+k2.*dz/2);
        k4       = fiber.Gamma.*(fiber.emi_sig.*N2-fiber.abs_sig.*(fiber.NT-N2) - fiber.Loss).*...
            (PSDdz1+k3.*dz);
        
        PSDdz2 = (PSDdz1 + dz/6*(k1 + 2*k2 + 2*k3 + k4)).*[ones(1,sys.nt),pumps.FWDchannels];
        pumps.Power(pumps.FWDchannels) = PSDdz2(logical([zeros(1,sys.nt),pumps.FWDchannels]));
        signal.PSD = PSDdz2(1:sys.nt); signal.Power = sum(signal.PSD);
    end
    fiber.G = (fiber.Gamma(1:sys.nt).*((fiber.emi_sig(1:sys.nt)+fiber.abs_sig(1:sys.nt)).*N2 - fiber.abs_sig(1:sys.nt).*fiber.NT) - fiber.Loss);
    %fiber.G = 1./(dz).*log((PSDdz2(1:sys.nt)+eps)./(PSDdz1(1:sys.nt)+eps)); % Calculate G using eps to combat machine error
    fiber.G(isnan(fiber.G)) = 00; % If gain is nan that means it was log(0/0) so set gain to 0
    
end

