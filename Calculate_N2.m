function N2 = Calculate_N2(sys,signal,pumps,ASE,fiber)
    run('units.m')
    PSDdz = [signal.PSD, pumps.PSD];
    lambda = [2*pi*c./(sys.w+2*pi*c./(sys.lambda0)),pumps.lambda];
    R12 = sum(fiber.Gamma.*fiber.abs_sig.*PSDdz./(hbar.*(2*pi*c./lambda).*fiber.aCore^2.*pi);
    R21 = sum(fiber.Gamma.*fiber.emi_sig.*PSDdz./(hbar.*(2*pi*c./lambda).*fiber.aCore^2.*pi);
    N2 = R12./(R12+R21+fiber.tau);
    
    R12 = sum(FWDpmp.Gamma.*FWDpmp.abs_sig.*FWDpmp.channels.*FWDpmp.PSD./(hbar.*(2*pi*c./FWDpmp.L).*fiber.Ac)) ...
        + sum(BKWpmp.Gamma.*BKWpmp.abs_sig.*BKWpmp.channels.*BKWpmp.PSD./(hbar.*(2*pi*c./BKWpmp.L).*fiber.Ac));
    R21 = sum(FWDpmp.Gamma.*FWDpmp.emi_sig.*FWDpmp.channels.*FWDpmp.PSD./(hbar.*(2*pi*c./FWDpmp.L).*fiber.Ac)) ...
        + sum(BKWpmp.Gamma.*BKWpmp.emi_sig.*BKWpmp.channels.*BKWpmp.PSD./(hbar.*(2*pi*c./BKWpmp.L).*fiber.Ac));
   
    
    W12 = sum(signal.Gamma.*fiber.abs_sig.*signal.channels.*signal.PSD./(hbar.*(signal.W).*fiber.Ac));
    W21 = sum(signal.Gamma.*fiber.emi_sig.*signal.channels.*signal.PSD./(hbar.*(signal.W).*fiber.Ac));
    
%     pumpFWDlam = 2*pi*c./(signal.w+pumpFWD.lambda0);
%     R12numer = 1./(hbar.*2*pi*c*fiber.Ac).*sum(pumpFWD.Gamma.*fiber.abs_sig.*pumpFWD.channels.*pumpFWD.PSD);
%     R12denom = 1./(hbar.*2*pi*c*fiber.Ac).*sum(pumpFWD.Gamma.*(fiber.abs_sig+fiber.emi_sig).*pumpFWD.channels
%     
    N2 = (R12 + W12)./(R12 + R21 + W12 + W21 + 1./fiber.tau).*fiber.NT;
    

end