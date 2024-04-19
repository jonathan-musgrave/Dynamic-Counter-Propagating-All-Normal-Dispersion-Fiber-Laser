function sys = plot_FiberPropResults(sys,signal,pumps,fiber)
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
            sys.fignum = sys.fignum+1;
   

end
