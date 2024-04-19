% Simulation of Bidrection fiber laser based on accurate modeling of the
% rate equations. 
% @author Jonathan Musgrave 2023
clear all
addpath([cd,'\Utils']);
addpath([cd,'\Plotting Function'])
run('units.m'); % Get units
sys = struct; passive_1 = struct; YDF = struct; passive_2 = struct; CWsignal = struct; CCWsignal = struct; pumps = struct;% initialize struct varibles
% Set system parameters
sys.PlotResults = 'no'; % We do not want to plot after every fiber propagation
sys.nt = 2^12;
sys.frep = 48*MHz; sys.nRT =699; sys.lambda0 =1075*nm;
step_size = 0.005;
sys.WinT = 30*ps; 
[sys, ~, ~, ~] = Fiber_Prop_Lz_v1(sys, CWsignal, passive_2, pumps); % Initialize

% Define the physical system parameters
passive_1.Lz = 1*m; passive_1.Type = 'passive';
YDF.Lz = 2.2*m; YDF.Type = 'gain'; YDF.RareEarthDopant = 'yb'; YDF.NT = 15e25; % Matching with Liekki Absorption
YDF.DoubleCladded = 'yes'; YDF.aCladding = 64*um;
passive_2.Lz = 1*m; passive_2.Type = 'passive';

% Define the CW signal and CCW signal 
CWsignal.Direction  = 'f'; % Direction is arbitrary for the signals but helps with organization with pumps
CCWsignal.Direction = 'b';
CWsignal.E0 = .25*pJ; CCWsignal.E0 = .25*pJ;
CWsignal.T0 = 200*fs; CCWsignal.T0 = 200*fs;

% Define the pumps
pumps.Direction = {'f','b'}; % 
pumps.lambda    = [976,976]*nm;
pumps.Power0    = [2.3,2.3]*W.*0 +[0.25,0.35]+[0.00,+.0];

% Define the waveplates
% Best result so far is 9*deg, -.3pi and .1 filter transmission

theta1= 83.15/180*pi;%9*deg;                   % Polarization Angle Of Analyzer(pi) 0.625
ph1=-.91*pi;

theta2=theta1;%-0.0003;                   % Polarization Angle Of Analyzer(pi) 0.625
ph2= ph1-0.00225;

[PC1CW,PC2CW,PC1CCW,PC2CCW] = get_Waveplates(theta1,theta2,ph1,ph2);

% Define the Spectral filter
BW_Filter = 2.62*2*pi*THz.*1.7 +0*2*pi*THz; % Given in frequency
BW_Filter_lam = 2*pi*c./(2*pi*c./sys.lambda0 - BW_Filter/2) - 2*pi*c./(2*pi*c./sys.lambda0 + BW_Filter/2);
Filt_trans = .4;
sys.SpectralFilter = Filt_trans.*exp(-2*log(2)*(sys.w./(BW_Filter)).^2);

% Calculate the total GDD


%% Begin Propagation Routine
CWuuOut_3d  = zeros(sys.nRT,sys.nt); 
CCWuuOut_3d = zeros(sys.nRT,sys.nt);
CWUUOut_3d = zeros(sys.nRT,sys.nt);
CCWUUOut_3d = zeros(sys.nRT,sys.nt);


dtShift_CW_RT = zeros(sys.nRT,1);
dtShift_CCW_RT = zeros(sys.nRT,1);
sys.OutputLength = 0; % for output
a = nan;

mod_PumpPower = 'yes';
Recenter = 'no';
pumpPowerMod = [300*mW,0]; % Amount of mod +- for CW and CCW directino
power0 = pumps.Power0;
nRTmod = 100; % number of roundtrips between power modulation
n = 1;
ncount = 1;
ncount1 = 1;
ncount2 = 1;
ncount3 = 1;
ncount4 = 1;
sys.GDD = (passive_1.Lz+passive_2.Lz+YDF.Lz).*20*fs.^2./mm;
for i = 1:sys.nRT % Do n Round Trip propagations
    sys.fignum = 1;
    % Propagate both through passive Fiber
    [sys, CWsignal, ~, pumps] = Fiber_Prop_Lz_v1(sys, CWsignal, passive_1, pumps, step_size*2);
    sys.fignum = sys.fignum+1;
    [sys, CCWsignal, ~, pumps] = Fiber_Prop_Lz_v1(sys, CCWsignal, passive_2, pumps, step_size*2); % Note the step size is twice here to speed up code since the NL phase is going to be low before gain section.
    sys.fignum = sys.fignum+1;
    
     if mod(i+5,nRTmod)==0 && (i+5)/nRTmod<3
        if ncount1 == 1
        ICWP1zNeg = CWsignal.Iw_3d; %Passive fiber 1 for CW
        ICCWP2zNeg = CCWsignal.Iw_3d; %Passive fiber 2 for CCW
        ncount1 = ncount1+1;
        else
        ICWP1zPos = CWsignal.Iw_3d; %Passive fiber 1 for CW
        ICCWP2zPos = CCWsignal.Iw_3d; %Passive fiber 2 for CCW
        end
    end
    
    % Propagate both through Gain Fiber
    % CWsignal.InterCavityPower = CCWsignal.Power; % Modify the saturation power by the internal energy of the counter propagating signal
    % CCWsignal.InterCavityPower = CWsignal.Power; % Modify the saturation power by the internal energy of the counter propagating signal
    [sys, CWsignal, ~, pumps] = Fiber_Prop_Lz_v1(sys, CWsignal, YDF, pumps, step_size);
    sys.fignum = sys.fignum+1;
    pumps.Direction = flip(pumps.Direction); % flip the diretion of the pumps for the CCW propagation
    [sys, CCWsignal, ~, pumps] = Fiber_Prop_Lz_v1(sys, CCWsignal, YDF, pumps, step_size);    
    sys.fignum = sys.fignum+1;
    pumps.Direction = flip(pumps.Direction); % flip the diretion of the pumps for the next round trip
    

    
    if mod(i+5,nRTmod)==0 && (i+5)/nRTmod<3
        if ncount == 1
        ICWGzNeg = CWsignal.Iw_3d;
        ICCWGzNeg = CCWsignal.Iw_3d;
        G_CWNeg = CWsignal.Gz_3d; % Sat Gain after YDF for CW
        G_CCWNeg = CCWsignal.Gz_3d; % Sat Gain after YDF for CCW
        ncount = ncount+1;
        else
        ICWGzPos = CWsignal.Iw_3d;
        ICCWGzPos = CCWsignal.Iw_3d;
        G_CWPos = CWsignal.Gz_3d; % Sat Gain after YDF for CW
        G_CCWPos = CCWsignal.Gz_3d; % Sat Gain after YDF for CCW
        end
    end
    
    % Propagate both through second passive fiber
    [sys, CWsignal, fiber_2, pumps] = Fiber_Prop_Lz_v1(sys, CWsignal, passive_2, pumps, step_size);
    G_RT2 = CWsignal.Gz_3d; % Returns 2d array same size as signal.uu_3d
    sys.fignum = sys.fignum+1;
    [sys, CCWsignal, fiber_1, pumps] = Fiber_Prop_Lz_v1(sys, CCWsignal, passive_1, pumps, step_size);
    G_RT1 = CCWsignal.Gz_3d; % Returns 2d array same size as signal.uu_3d
    sys.fignum = sys.fignum+1;

    if mod(i+5,nRTmod)==0 && (i+5)/nRTmod<3
        if ncount2 == 1
        ICWP2zNeg = CWsignal.Iw_3d; %Passive fiber 2 for CW
        ICCWP1zNeg = CCWsignal.Iw_3d; %Passive fiber 1 for CCW
        ncount2 = ncount2+1;
        else
        ICWP2zPos = CWsignal.Iw_3d; %Passive fiber 2 for CW
        ICCWP1zPos = CCWsignal.Iw_3d; %Passive fiber 1 for CCW
        end
    end

    
    % now for Waveplate Propagation and output coupling of h
    % CW signal;
    [CWsignal.uu] = PC1CW*[CWsignal.h;CWsignal.v]; % Apply Waveplates
    if strcmpi(Recenter,'yes');
        CWsignal.uu = fftshift(circshift(CWsignal.uu,[0,-SHIFTINX]),2);
    end
    CWsignal.uuOut = CWsignal.uu(1,:); % Get Current Output of laser H pol
    CWsignal.V = sys.SpectralFilter.*IFFT(CWsignal.uu(2,:),sys.nt,sys.dt,sys.frep);  % Apply spectral filter
    [CWsignal.uu] =  PC2CW*[FFT(CWsignal.V,sys.nt,sys.dt,sys.frep);zeros(1,sys.nt)];  % Apply secon set of Waveplates
    CWsignal.h = CWsignal.uu(1,:); CWsignal.v = CWsignal.uu(2,:); % Update CW signal; (H and V are updated automattically in Fiber_Prop_Lz)
    
    if mod(i+5,nRTmod)==0 && (i+5)/nRTmod<3
        if ncount3 == 1
        ICWAOCzNeg = CWsignal.V; %Passive fiber 2 for CW
        ncount3 = ncount3+1;
        else
        ICWAOCzPos = CWsignal.V; %Passive fiber 2 for CW
        end
    end    
    
    [CCWsignal.uu] = PC1CCW*[CCWsignal.h;CCWsignal.v]; % Apply Waveplates
    if strcmpi(Recenter,'yes');
        CCWsignal.uu = fftshift(circshift(CCWsignal.uu,[0,-SHIFTINX]),2);
    end
    CCWsignal.uuOut = CCWsignal.uu(1,:); % Get Current Output of laser H pol
    CCWsignal.V = sys.SpectralFilter.*IFFT(CCWsignal.uu(2,:),sys.nt,sys.dt,sys.frep);  % Apply spectral filter
    [CCWsignal.uu] =  PC2CCW*[FFT(CCWsignal.V,sys.nt,sys.dt,sys.frep);zeros(1,sys.nt)];  % Apply secon set of Waveplates
    CCWsignal.h = CCWsignal.uu(1,:); CCWsignal.v = CCWsignal.uu(2,:); % Update CW signal; (H and V are updated automattically in Fiber_Prop_Lz)
    
     if mod(i+5,nRTmod)==0 && (i+5)/nRTmod<3
        if ncount4 == 1
        ICCWAOCzNeg = CCWsignal.V; %Passive fiber 2 for CW
        ncount4 = ncount4+1;
        else
        ICCWAOCzPos = CCWsignal.V; %Passive fiber 2 for CW
        end
    end   
   
    %% Round trip count
    fprintf(['Round Trip \n', num2str(i)]);
    
    CWuuOut_3d(i,:) = CWsignal.uuOut; 
    CCWuuOut_3d(i,:) = CCWsignal.uuOut;
    CWUUOut_3d(i,:) = IFFT(CWuuOut_3d(i,:),sys.nt,sys.dt,sys.frep);
    CCWUUOut_3d(i,:) = IFFT(CCWuuOut_3d(i,:),sys.nt,sys.dt,sys.frep);
    IwCW  = CWUUOut_3d(i,:).*conj(CWUUOut_3d(i,:));
    IwCCW  = CCWUUOut_3d(i,:).*conj(CCWUUOut_3d(i,:));
    ItCW  = CWuuOut_3d(i,:).*conj(CWuuOut_3d(i,:));
    ItCCW = CCWuuOut_3d(i,:).*conj(CCWuuOut_3d(i,:));
    CWEng(i) = sum(IwCW)./sys.frep; 
    CCWEng(i) =   sum(IwCCW)./sys.frep;
    
    COM_tCW(i) =calc_COM(sys.t,ItCW);
    COM_tCCW(i) =calc_COM(sys.t,ItCCW);
    if i>2 % dCOM_t/dRT = delta_frep from sys.frep
        % (d/dn(COM_t(n)) +1/sys.frep)^(-1) = frep
        % Calculate slope
        y1CW = COM_tCW(i-1);
        y2CW = COM_tCW(i);
        dCOM_dn = y2CW - y1CW;
        frepCW(i) = 1./(dCOM_dn/2/pi +1./sys.frep);
        
        
        y1CCW = COM_tCCW(i-1);
        y2CCW = COM_tCCW(i);
        dCOM_dn = y2CCW - y1CCW;
        % Recenter if pulse is walking off too much
        if y2CW > sys.WinT./2 || y2CCW>sys.WinT./2
            Recenter = 'yes';
            [~,SHIFTINX] = min(abs(sys.t - max([y2CW,y2CCW])));
        else
            Recenter = 'no';
        end
        frepCCW(i) = 1./(dCOM_dn/2/pi +1./sys.frep);
        dFrep(i) = frepCW(i) -frepCCW(i);
        LW = 2;
        figure(sys.fignum-1);clf;
        subplot(2,1,1);
        plot(COM_tCW(1:i)./fs,'-r','linewidth',LW)
        hold on;
        plot(COM_tCCW(1:i)./fs,'-b','linewidth',LW);
        xlabel('RT')
        ylabel('t_COM(n)')
        
        subplot(2,1,2);
        hold on;
        plot(frepCW(1:i)./MHz,'-r','linewidth',LW)
        plot(frepCW(1:i)./MHz,'or')
        plot(frepCCW(1:i)./MHz,'-b','linewidth',LW);
        plot(frepCCW(1:i)./MHz,'ob');
        xlabel('RT')
        ylabel('frep(n)')
        yyaxis('right');
        plot(dFrep(1:i));
        ylabel('\Delta f_{rep}')
        set(gca,'YColor',[0,1,0]./2)
        ylim([-5 5])
        %ylim(dFrep(i)+[-5,5]);
        drawnow
    end
    
    
        
        
        
    if (strcmpi(mod_PumpPower,'yes') && mod(i+1,nRTmod)==0)
       pumps.Power0 = power0 + n*pumpPowerMod;
       %n = n*-1; %Switch mod direction
       n = abs(n-1); %Switch mod direction
    end
    
    saveh=figure(sys.fignum);clf; hold on;
    subplot(2,3,1);hold on;
    plot(1:i,CWEng(1:i)./nJ,'-r')
    plot(1:i,CWEng(1:i)./nJ,'or')
    plot(1:i,CCWEng(1:i)./nJ,'-b')
    plot(1:i,CCWEng(1:i)./nJ,'ob')
    ylabel('Pulse Energy (nJ)')
    title({['P_{CW}: ', num2str(pumps.Power0(1))] ['P_{CCW}: ', num2str(pumps.Power0(2))]})
    
    subplot(2,3,2);hold on
    W2dBm = @(x) x; % For linear plotting
    
    y1 = W2dBm(fftshift(IwCW));
    y2 = W2dBm(fftshift(IwCCW));
    plot(2*pi*c./(fftshift(sys.w)+2*pi*c./sys.lambda0)./nm,y1 - max(y1),'-r')
    plot(2*pi*c./(fftshift(sys.w)+2*pi*c./sys.lambda0)./nm,y2 - max(y2),'--b')
    % ylim([-50,0])
    plot(2*pi*c./(fftshift(sys.w)+2*pi*c./sys.lambda0)./nm,y1./max(y1),'-r')
    plot(2*pi*c./(fftshift(sys.w)+2*pi*c./sys.lambda0)./nm,y2./max(y2),'--b')
    ylim([0,1.2]);
    xlim([1050,1100])
    ylabel('Pulse Energy (nJ)')
    xlabel('\lambda (nm)')
    drawnow;
    
    subplot(2,3,3);hold on
    y1 = ItCW./max(ItCW);
    y2 = ItCCW./max(ItCCW);
    plot(sys.t./ps,y1,'-r')
    plot(sys.t./ps,y2,'--b')
    ylim([0,1.2])
    ylabel('Pulse Energy (nJ)')
    xlabel('t (ps)')
    
    
    % Calculate curent delta frep
    lamCOM_CW(i) = calc_CarrierFreq(fftshift((IwCW).^(1/2))',2*pi*c./(fftshift(sys.w)+2*pi*c./sys.lambda0));
    lamCOM_CCW(i) = calc_CarrierFreq(fftshift((IwCCW).^(1/2))',2*pi*c./(fftshift(sys.w)+2*pi*c./sys.lambda0));
    
    del_CWw =  -(2*pi*c./(sys.lambda0) - 2*pi*c./(lamCOM_CW(i))); % red in fig
    del_CCWw =  -(2*pi*c./(sys.lambda0) - 2*pi*c./(lamCOM_CCW(i))); % blue in fig
    del_centerw = (del_CWw - del_CCWw);
    dDfrep(i) = -sys.frep.^2*(sys.GDD.*del_centerw/2/pi); %changed this term
    
    subplot(2,3,[4:6]);hold on
    y1 = lamCOM_CW(1:i)./nm;
    y2 = lamCOM_CCW(1:i)./nm;
    plot(1:i,y1,'-r');
    plot(1:i,y1,'or')
    plot(1:i,y2,'-b');
    plot(1:i,y2,'ob');
    ylabel('\lambda_{COM} (nm)')
    
    yyaxis('right'); set(gca,'YCOLOR',[0,0.5,0]);
    plot(1:i,(dDfrep(1:i)),'Color',[0,0.5,0])
    plot(1:i,(dDfrep(1:i)),'o','Color',[0,0.5,0])
    ylim([-5 5])

 
%     ylim(dDfrep_2(i) + [-30,30])
    ylabel('\Delta f_{rep}')
    drawnow;
    
        %y2
   if mod(i+3,nRTmod)==0 && (i+3)/nRTmod<4
        
         if i<nRTmod
            annstr = sprintf('Before Modulation'); % annotation text
            annpos = [0.5 0.89 0.1 0.1]; % annotation position in figure coordinates
            ha = annotation('textbox',annpos,'string',annstr,'FontSize',14);
            ha.HorizontalAlignment = 'center';            
         elseif n<0
            annstr = sprintf('Positive Pump Power Modulation'); % annotation text
            annpos = [0.6 0.89 0.1 0.1]; % annotation position in figure coordinates
            ha = annotation('textbox',annpos,'string',annstr,'FontSize',14);
            ha.HorizontalAlignment = 'center';  
         elseif n>0
            annstr = sprintf('Negative Pump Power Modulation'); % annotation text
            annpos = [0.6 0.89 0.1 0.1]; % annotation position in figure coordinates
            ha = annotation('textbox',annpos,'string',annstr,'FontSize',14);
            ha.HorizontalAlignment = 'center';
         end
          saveas(saveh,sprintf('FIG%d.png',(i+3)/nRTmod)); % will create FIG1, FIG2,...
    end
       
    
end







function [PC1CW,PC2CW,PC1CCW,PC2CCW] = get_Waveplates(theta1,theta2,ph1,ph2)
    PC1CW=[cos(theta1),-sin(theta1);sin(theta1),cos(theta1)]*[exp(-1i*ph1/2),0;0,exp(1i*ph1/2)]*[cos(theta1),sin(theta1);-sin(theta1),cos(theta1)];
    PC2CW=[cos(theta2),-sin(theta2);sin(theta2),cos(theta2)]*[exp(-1i*ph2/2),0;0,exp(1i*ph2/2)]*[cos(theta2),sin(theta2);-sin(theta2),cos(theta2)];
    PC1CCW=[cos(-theta2),-sin(-theta2);sin(-theta2),cos(-theta2)]*[exp(1i*ph2/2),0;0,exp(-1i*ph2/2)]*[cos(-theta2),sin(-theta2);-sin(-theta2),cos(-theta2)];
    PC2CCW=[cos(-theta1),-sin(-theta1);sin(-theta1),cos(-theta1)]*[exp(1i*ph1/2),0;0,exp(-1i*ph1/2)]*[cos(-theta1),sin(-theta1);-sin(-theta1),cos(-theta1)];
end