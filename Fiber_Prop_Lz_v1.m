function [sys, signal, fiber, pumps] = Fiber_Prop_Lz_v1(sys, signal, fiber, pumps, step_size)
    addpath([cd,'\Utils'])
    run('units.m')
    % Check for sys fields and set default if not there
    if ~exist('sys','var'); sys = struct; end
    if ~isfield(sys,'nt'), sys.nt = 2^10; end % Default sampling
    if ~isfield(sys,'WinT'), sys.WinT = 15*ps; end % Default one sided time window
        sys.dt = 2*sys.WinT./sys.nt;
        sys.dw = 2*pi./(sys.nt*sys.dt); 
    if ~isfield(sys,'lambda0'), sys.lambda0 = 1030*nm; end % Default sampling
    if ~isfield(sys,'t'),    sys.t = (-sys.nt/2:(sys.nt/2-1))*sys.dt; end 
    if ~isfield(sys,'w'),    sys.w = fftshift(2*pi/(sys.nt*sys.dt)*(-sys.nt/2:sys.nt/2-1)); end
    if ~isfield(sys,'frep'), sys.frep = 80*MHz; end % Default repetition rate 
    if ~isfield(sys,'PlotResults'), sys.PlotResults = 'yes';end % Plot output of system including 3d propagation
    if ~isfield(sys,'fignum'), sys.fignum = 1; end % Define the fignum start
    if ~isfield(sys,'BoundaryCondition'), sys.BoundaryCondition = ones(1,sys.nt);end
    if ~isfield(sys,'MaxIter'), sys.MaxIter = 20;end
    if ~isfield(sys,'tol'),sys.tol = .0001; end % RMSE Tolerance of shooting method
        
    % Check for signal fields and set to default if not there
    if ~exist('signal','var'); signal = struct; end
    if ~isfield(signal,'T0'), signal.T0 = .3*ps; end % Default FWHM of pulse intensity (in time)
    if ~isfield(signal,'E0'), signal.E0 = .25*nJ; end % Default pulse energy 
    if ~isfield(signal,'Power0'), signal.Power0 = signal.E0.*sys.frep; end
    if ~isfield(signal,'uu0'), signal.uu0 = exp(-4*log(2)*(sys.t./signal.T0).^2); 
        signal.uu0 = signal.uu0./(sqrt(sum(abs(signal.uu0).^2.*sys.dt./signal.Power0.*sys.frep))); % Set the power
        signal.UU0 = IFFT(signal.uu0,sys.nt,sys.dt,sys.frep); 
        signal.PSD0 = abs(signal.UU0).^2;
        signal.uu = signal.uu0; signal.UU = signal.UU0;
    end % Default is gaussian signal
    if ~isfield(signal,'h')||~isfield(signal,'v') % Default polarization is 1/sqrt(2)*[uu0;uu0] i.e P45
        signal.h0 = signal.uu0./sqrt(2); signal.h = signal.h0;
        signal.v0 = signal.uu0./sqrt(2).*exp(1j*2*pi*0); signal.v = signal.v0;
    end
    if ~isfield(signal,'Direction'), signal.Direction = {'f'}; end
    if ~isfield(signal,'InterCavityPower'), signal.InterCavityPower = 0; end % Modify the intercavity Power usually from counter propagating signals
    % Multiple signals up to two can be defined for propagation
    signal.H =  IFFT(signal.h,sys.nt,sys.dt,sys.frep); signal.H0 = signal.H;
    signal.V =  IFFT(signal.v,sys.nt,sys.dt,sys.frep); signal.V0 = signal.V;
    signal.h0 = signal.h; 
    signal.v0 = signal.v;
    signal.PSD = abs(signal.H).^2 + abs(signal.V).^2;
    signal.Power = sum(signal.PSD);

        
    
    
    % Check for fiber fields and set default if not there
    if ~exist('fiber','var'); fiber = struct; end
    if ~isfield(fiber,'Raman'), fiber.Raman = 'off'; end % Default to no Raman
    if ~isfield(fiber,'Type'), fiber.Type = 'passive'; end % Default to passive fiber propagation
    if ~isfield(fiber,'betas'), fiber.betas = [23*fs^2./mm, 30*fs^3./mm]; end % Default Dispersion is HI-1060
    if ~isfield(fiber,'n2'), fiber.n2 = 2.19*10^(-20)*m^2./W; end
    if ~isfield(fiber,'NACore'),fiber.NACore = 0.13; end % Core NA
    if ~isfield(fiber,'aCore'), fiber.aCore = 5*um/2; end % Core Radius
    if ~isfield(fiber,'Lz'), fiber.Lz = 2.5*m; end % Fiber length
    if ~isfield(fiber,'alpha'), fiber.alpha = 0;end % fiber Loss in db/m (positive for loss)
    % fiber defaults for gain fiber
    if strcmpi(fiber.Type,'gain')
        if (~isfield(fiber,'NT')), fiber.NT = 1e25; end
        if ~isfield(fiber,'RareEarthDopant'),fiber.RareEarthDopant = 'yb';end
        if ~isfield(fiber,'DoubleCladded'),fiber.DoubleCladded = 'no';end
        if ~isfield(fiber,'tau'),fiber.tau = 1.4*ms; end % Relaxation rate
        if strcmpi(fiber.DoubleCladded,'yes') % If double cladded define parameter
            if ~isfield(fiber,'CladdingNA'),fiber.CladdingNA = 0.2; end
            if ~isfield(fiber,'aCladding'), fiber.aCladding  = 125*um/2; end
        end
    
        if ~exist('pumps','var') && strcmpi(fiber.Type,'gain'), pumps = struct;end % Check if pumps are defined if not save as empty
        if ~isfield(pumps,'lambda'), pumps.lambda = [976*nm]; end
        if ~isfield(pumps,'Power0'), pumps.Power0 = [0.6]*W; end
        if ~isfield(pumps,'Direction'), pumps.Direction = {'f'}; end
        if ~isfield(pumps,'np'), pumps.np = length(pumps.lambda); end
        pumps.FWDchannels = strcmpi(pumps.Direction,'f');
        pumps.BKWchannels = strcmpi(pumps.Direction,'b');
        pumps.Power = pumps.Power0;
    end
    if ~exist('pumps','var') && strcmpi(fiber.Type,'passive'), pumps = [];end % If passive make empty
    
    if ~exist('step_size','var') % Default Step size
        step_size = 1*cm;
    end
    %% Calculate loss profile
    fiber.Loss = 1-0.^(-fiber.alpha./10); % Convert db/m to percentage/m
    %% Define the dispersion profile operator
    fiber.dispersion = 0;
    for i = 1:length(fiber.betas) % When defining fiber the dispersion operator is saved as an array
        % Calculated with an fftshifted omega
       fiber.dispersion =  (1./(factorial(i+1)).*fiber.betas(i).*...
           sys.w.^(i+1))+fiber.dispersion; % Taylor expansion definition without the dz added
    end
    %% Define the effctive area and the resulting nonlinearity parameter
    if strcmpi(fiber.Type,'passive') 
        V = [sys.w+2*pi*c./sys.lambda0]./c.*fiber.aCore.*fiber.NACore; % V num of [signal,pumps]
    else % if not passive calculate the V number of the pumps as well
        V = [sys.w+2*pi*c./sys.lambda0,2*pi*c./(pumps.lambda)]./c.*fiber.aCore.*fiber.NACore; % V num of [signal,pumps]
    end  
    Wi = fiber.aCore*(0.616 + 1.66./V.^(1.5) + 0.987./V.^6);
    fiber.Aeff  = pi*Wi.^2; % gaussian mode area approximation
    fiber.gamma = 2*pi./sys.lambda0.*fiber.n2./fiber.Aeff(1); % Just calculates the center wavelength nonlinearirt
    %% Calculate the mode overlap with the core
    fiber.Gamma = 1-exp(-2*(fiber.aCore./Wi).^2);
    if strcmpi(fiber.Type,'gain') % If gain fiber check if double cladded
        if strcmpi(fiber.DoubleCladded,'yes') % If double cladded the pump mode overlap is modulated (signal still travels in core)
           fiber.Gamma(end-pumps.np+1:end) = fiber.Gamma(end-pumps.np+1:end).*(fiber.aCore./fiber.aCladding).^2; % modulated by ratio of areas 
        end
    end
    
    %% Define Raman term if necessary
    if strcmpi(fiber.Raman,'on')
        % Raman Response Updated (not the same as paper but is the updated versions  from Agrawal et. al.) 
        fiber.fR = 0.235;  fb=0.21;
        tau1 = 12.2.*fs; tau2 = 32.*fs; tau_b = 96.*fs; 
        fiber.h = (1-fb)*(tau1.^2 +tau2.^2)./(tau1*tau2^2)*exp(-sys.t/tau2).*sin(sys.t/tau1) ...
            + fb*((2*tau_b-sys.t)/tau_b^2.*exp(-sys.t/tau_b));
        fiber.h(sys.t<0) = 0; % Causality (Heavy side step)
    end
    %%%%%%%%%% IF FIBER IS GAIN DEFINE THE RESULTING PARAMETERS %%%%%%%%%%%
    if strcmpi(fiber.Type,'gain')
    %% Define the absorption and emission cross section if 
        switch fiber.RareEarthDopant
            case 'yb'
                % Emission and absorption cross sections for Yb Doped fiber for
                % the signal wavelength
                [fiber.emi_sig, fiber.abs_sig] = getYDFSpectroData_v2([(sys.w+2*pi*c./sys.lambda0),2*pi*c./pumps.lambda]); % Get the Spectroscopy data interpolated over the frequency axis
                
        end
        
        %% Calculate the saturation energy
        fiber.Esat = fiber.Aeff(1:end-pumps.np).*hbar.*([sys.w+2*pi*c./sys.lambda0])./(2*pi.*(fiber.emi_sig(1:end-pumps.np)+fiber.abs_sig(1:end-pumps.np)));
        fiber.Esat(fiber.Esat > 100*uJ) = max(fiber.Esat(fiber.Esat<100*uJ)); % Editing so there are no inf values  
        fiber.Psat = fiber.Esat(1).*sys.frep; % Saturation power 
    end
    
    
    %% MAIN PROPAGATION ROUTINES
    zi = 0; % current z
    n = 1; % Current save
    nz = floor(fiber.Lz./step_size)+1; % total nz steps +1 for 0(no Adaptive step as of 9/28/2023)
    signal.uz_3d = zeros(nz,sys.nt); 
    signal.Uz_3d = zeros(nz,sys.nt); 
    sys.z = zeros(nz,1); n = n+1; % increase n. z(1) = 0;
    % Choose propagation routine based on the fiber type (defined when
    % fiber is described)
    switch fiber.Type
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case 'passive' % Passive propagation
            % prealocate
            signal.hz_3d = zeros(nz,sys.nt); signal.vz_3d = zeros(nz,sys.nt);
            signal.hz_3d(1,:) = signal.h; signal.vz_3d(1,:) = signal.v;
            signal.Hz_3d = zeros(nz,sys.nt); signal.Vz_3d = zeros(nz,sys.nt);
            signal.Hz_3d(1,:) = signal.H; signal.Vz_3d(1,:) = signal.V;
            signal.Powerz = zeros(nz,1);  signal.Powerz = sum(signal.PSD);
            while zi < fiber.Lz % Using while if adaptive step is  to be used as well
                dz = step_size;
                Gsat = 0; 
                signal = PropagateCNLSE_dz(sys,signal, fiber, Gsat,  dz);
                signal.hz_3d(n,:) = signal.h; signal.vz_3d(n,:) = signal.v;
                signal.Hz_3d(n,:) = signal.H; signal.Vz_3d(n,:) = signal.V;
                signal.Powerz(n,:) = sum(signal.PSD);
                zi = zi+dz;
                sys.z(n) = zi;
                n = n+1;
            end
            signal.Ez = signal.Powerz./sys.frep; % Calculate pulse Energy
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case 'gain' % Gain fiber propagation using IVP shooting method. 
            % Preallocate variables
            signal.hz_3d = zeros(nz,sys.nt); signal.vz_3d = zeros(nz,sys.nt);
            signal.hz_3d(1,:) = signal.h; signal.vz_3d(1,:) = signal.v;
            signal.Hz_3d = zeros(nz,sys.nt); signal.Vz_3d = zeros(nz,sys.nt);
            signal.Hz_3d(1,:) = signal.H; signal.Vz_3d(1,:) = signal.V;
            signal.Powerz = zeros(nz,1);  signal.Powerz = sum(signal.PSD);
            fiber.N2z = zeros(nz-1,1); % minus 1 because its N2z/m;
            
            % Get forward and backward pumps
            pumps.Powerz = zeros(nz,pumps.np);
            if max(pumps.FWDchannels), pumps.Powerz(1,pumps.FWDchannels) = pumps.Power0(pumps.FWDchannels);end % Get forward pump Powers
            if max(pumps.BKWchannels), pumps.Powerz(:,pumps.BKWchannels) = repmat(pumps.Power0(pumps.BKWchannels),[nz,1]);end % Get backward pump Power initial guess for IVP 
            
            
            if ~max(pumps.BKWchannels); % if no Backwardward channels do forward prop
                % Begin forward propagation
                while zi < fiber.Lz 
                    dz = step_size;
                    % Solve the Rate equation for dz (solves N2)
                    [sys, signal, pumps, fiber] = Rate_Equation_NoASEdz(sys, signal, pumps, fiber, dz);
                    % Save pump power
                    pumps.Powerz(n,:) = pumps.Power;
                    signal.Powerz(n) = signal.Power;
                    Gsat = fiber.G.*exp(-((signal.Power+signal.InterCavityPower)./fiber.Psat));

                    signal.Gz_3d((n-1),:) = Gsat;
                    fiber.N2z(n-1) = fiber.N2;
                    % Calculate G
                    signal = PropagateCNLSE_dz(sys, signal, fiber, Gsat,  dz);
                    signal.hz_3d(n,:) = signal.h; signal.vz_3d(n,:) = signal.v;
                    signal.Hz_3d(n,:) = signal.H; signal.Vz_3d(n,:) = signal.V;

                    zi = zi+dz;
                    sys.z(n) = zi;
                    n = n+1;
                end
                
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            else % Iterative propagation for IVP problem
                k = 1;% iter value
                OutputLength = 0; % For iteration display
                % Do a Backwards Propagation to solve the power
                % distribution of the backwards pump
                
                while k <= sys.MaxIter
                    n = 1;
                    signal.h = signal.h0; signal.v = signal.v0;
                    signal.hz_3d(n,:) = signal.h; signal.vz_3d(n,:) = signal.v;
                    signal.Hz_3d(n,:) = signal.H0; signal.Vz_3d(n,:) = signal.V0;
                    signal.PSD = signal.PSD0; signal.Powerz(n) = sum(signal.PSD);
                    signal.Power = signal.Power0;
                    
                    
                    n = n+1;
                    zi = 0;
                    % Do a forward Propagation
                    while zi < (fiber.Lz - step_size/2)
                        dz = step_size;
                        pumps.Power(pumps.BKWchannels) = pumps.Powerz(n,pumps.BKWchannels); % Get the BKW pump solution from backward propagation (or initial solution guess)
                        
                        % Solve the Rate equation for dz (solves N2)
                        [sys, signal, pumps, fiber] = Rate_Equation_NoASEdz(sys, signal, pumps, fiber, dz);
                        % Save pump power
                        pumps.Powerz(n,pumps.FWDchannels) = pumps.Power(pumps.FWDchannels);
                        signal.Powerz(n) = signal.Power;
                        fiber.G3d(n,:) = fiber.G;
                        Gsat = fiber.G.*exp(-((signal.Power+signal.InterCavityPower)./fiber.Psat));

                        signal.Gz_3d((n-1),:) = Gsat;
                        fiber.N2z(n-1) = fiber.N2;
                        % Calculate G
                        signal = PropagateCNLSE_dz(sys, signal, fiber, Gsat,  dz);
                        signal.hz_3d(n,:) = signal.h; signal.vz_3d(n,:) = signal.v;
                        signal.Hz_3d(n,:) = signal.H; signal.Vz_3d(n,:) = signal.V;
                        signal.Powerz(n) = signal.Power;
                            
                        zi = zi+dz;
                        sys.z(n) = zi;
                        n = n+1;
                    end
                    signalFWD_PSDz = abs(signal.Hz_3d).^2+abs(signal.Vz_3d).^2;
                    % plot_Results(sys,signal,pumps,fiber)
                    PSD_1 = signal.PSD; % Save output power of signal
                    if k>1
                        eps1 = PSD_2 - PSD_1;
                        RMSE = sqrt(sum(eps1.^2)./(sys.nt)); % Calculate RMSE error of the FWD vs Backward Prop
                        fprintf(repmat(['\b'],[1,OutputLength]))
                        Out1 = fprintf(['Iteration: ',num2str(k),'\n RMSE: ']);
                        Out2 = fprintf('%.4f',[RMSE]); 
                        
                        OutputLength = Out1+Out2 + fprintf('\n');
                        
                        if RMSE<sys.tol, sys.OutputLength = OutputLength; break; end
                    end
                    PSD_2 = PSD_1;
                    
                    % Do a Backwards Propagation to solve the power
                    % distribution of the backwards pump
                    n = n-1;
                    while zi > step_size/2
                        dz = -step_size;
                        % Solve the Rate equation for -dz (solves N2)
                        % Get forward solved values for N2, pump power and
                        % signal pSD
                        pumps.Power(pumps.FWDchannels) = pumps.Powerz(n-1,pumps.FWDchannels); % Get the FWD pump solution from forward propagation
                        signal.PSD = signalFWD_PSDz(n,:);
                        fiber.N2 = fiber.N2z(n-1);
                        [sys, signal, pumps, fiber] = Rate_Equation_NoASEdz(sys, signal, pumps, fiber, dz);
                        % Save pump power
                        pumps.Powerz(n-1,pumps.BKWchannels) = pumps.Power(pumps.BKWchannels); % update the backwards channels
                        
                        zi = zi+dz;
                        sys.z(n-1) = zi;
                        n = n-1; 
                    end
                    k = k+1;
                    
                    
                end
                    
            end
            
            
    end
    % Get the intensity in time and frequency
    signal.h = signal.hz_3d(end,:); signal.v = signal.vz_3d(end,:);
    signal.H = signal.Hz_3d(end,:); signal.V = signal.Vz_3d(end,:);
    signal.It_3d = abs(signal.hz_3d).^2+abs(signal.vz_3d).^2;
    signal.Iw_3d = abs(signal.Hz_3d).^2+abs(signal.Vz_3d).^2;
    signal.It = abs(signal.hz_3d(end,:)).^2+abs(signal.vz_3d(end,:)).^2;
    signal.Iw = abs(signal.Hz_3d(end,:)).^2+abs(signal.Vz_3d(end,:)).^2;
    
    if strcmpi(sys.PlotResults,'yes')
         plot_Results(sys,signal,pumps,fiber)
    end
        
end

function plot_Results(sys,signal,pumps,fiber)
        run('units.m')
        signal.It = ((abs(signal.hz_3d).^2+abs(signal.vz_3d).^2));
        signal.Iw = ((abs(signal.Hz_3d).^2+abs(signal.Vz_3d).^2));
        Fontname = 'Helvetica';
        LW = 3; FS = 10; % Linewidth and Fontsize
        zUnits = cm; PUnits = W; EUnits = nJ; tUnits = ps; wUnits = 2*pi*THz;
        zStr   = 'z (cm)'; PStr = 'Power (W)'; EStr = 'Energy (nJ)'; tStr = 'time (ps)';  wStr = 'THz';
        fig = figure(sys.fignum); clf; 
        fig.Position = [0       0         900         850];
        
        wCLIM = 10*log10(1e3*[min(min(signal.Iw(signal.Iw~=0))),max(max(signal.Iw))]); % conversion into dBm 
            wCLIM(1) = wCLIM(2) - 60;
            
        tCLIM = [0,max(max(signal.It))];
        
        subplot(5,2,1); hold on;
            plot(sys.t./tUnits,signal.It(end,:),'-r','linewidth',LW);
            xlabel(tStr); ylabel(tStr); ylim([0,max(signal.It(end,:)).*1.2]);
        subplot(5,2,2);hold on;
            plot(2*pi*c./(fftshift(sys.w)+2*pi*c./sys.lambda0)./nm,fftshift(signal.Iw(end,:))./max(signal.Iw(end,:)),'-r','linewidth',LW);
            xlabel('\lambda (nm)'); ylim([0,1.2]); ylabel('dBm')
            xlim([-20,20]+sys.lambda0./nm)
        subplot(5,2,[3,4]); hold on;grid on;
            plot(sys.z./zUnits, signal.Powerz./PUnits,'-r','linewidth',LW);
            plot(sys.z(1:20:end)./zUnits, signal.Powerz(1:20:end)./PUnits,'or','linewidth',LW); 
            ylabel(PStr);ylim([0,max(signal.Powerz)].*1.2);
            if strcmpi(fiber.Type,'gain')
                plot(sys.z./zUnits, (pumps.Powerz)./PUnits,'-b','linewidth',LW/2);
            	ylim([0,max(max(signal.Powerz),max(max(sum(pumps.Powerz,2))))].*1.2);
                yyaxis('right'); plot(sys.z(2:end)./zUnits,fiber.N2z,'--k'); ylabel('N2 (%)'); set(gca,'Ycolor','k');
                ylim([0,1])
            end
            xlim([0,fiber.Lz]./zUnits)
            
            
        % S-Polarized light plotting (i.e. H)
        s1 = subplot(5,2,5); hold on;
            imagesc(sys.z./zUnits,sys.t./tUnits, abs(signal.hz_3d').^2)
            xlabel(zStr); ylabel(tStr);
            xlim([0,fiber.Lz]./zUnits);
            ylim([min(sys.t),max(sys.t)]./tUnits);
            set(gca,'Clim',tCLIM)
        s2 = subplot(5,2,6); hold on;
            imagesc(sys.z./zUnits,fftshift(sys.w)./wUnits, 10*log10(1e3*fftshift(abs(signal.Hz_3d').^2,1)))
            xlabel(zStr); ylabel(wStr);
            xlim([0,fiber.Lz]./zUnits);
            ylim([min(sys.w),max(sys.w)]./wUnits);
            set(gca,'Clim',wCLIM)
        
        % P-Polarized light plotting (i.e. H)
        s3 = subplot(5,2,7); hold on;
            imagesc(sys.z./zUnits,sys.t./tUnits, abs(signal.vz_3d').^2)
            xlabel(zStr); ylabel(tStr);
            xlim([0,fiber.Lz]./zUnits);
            ylim([min(sys.t),max(sys.t)]./tUnits);
            set(gca,'Clim',tCLIM)
        s4 = subplot(5,2,8); hold on;
            imagesc(sys.z./zUnits,fftshift(sys.w)./wUnits, 10*log10(1e3*fftshift(abs(signal.Vz_3d').^2,1)))
            xlabel(zStr); ylabel(wStr);
            xlim([0,fiber.Lz]./zUnits);
            ylim([min(sys.w),max(sys.w)]./wUnits);
            set(gca,'Clim',wCLIM)
        s5 = subplot(5,2,9); hold on;
            imagesc(sys.z./zUnits,sys.t./tUnits, signal.It')
            xlabel(zStr); ylabel(tStr);
            xlim([0,fiber.Lz]./zUnits);
            ylim([min(sys.t),max(sys.t)]./tUnits);
            set(gca,'Clim',tCLIM);ct = colorbar;
        s6 = subplot(5,2,10); hold on;
            imagesc(sys.z./zUnits,fftshift(sys.w)./wUnits, 10*log10(1e3*fftshift(signal.Iw',1)))
            xlabel(zStr); ylabel(wStr);
            xlim([0,fiber.Lz]./zUnits);
            ylim([min(sys.w),max(sys.w)]./wUnits);
            set(gca,'Clim',wCLIM); cw = colorbar;
            colormap(s1,'hot');colormap(s3,'hot');colormap(s5,'hot');
            colormap(s2,'jet');colormap(s4,'jet');colormap(s6,'jet');
            
            ct.Position = [0.061770202020202,0.156145287631589,0.013785353535354,0.36856059472135];
            cw.Position =[0.934629629281936,0.137647058823527,0.012037037384731,0.369411762407102];
   

end
