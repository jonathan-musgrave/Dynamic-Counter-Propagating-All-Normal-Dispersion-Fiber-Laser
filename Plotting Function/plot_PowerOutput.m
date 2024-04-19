function sys = plot_PowerOutput(sys,signal,pumps,fiber,options)
if ~exist('options','var'), options = struct;end
if ~isfield(options,'units'), options.units = 'linear';end
run('units.m');
if strcmpi(options.units,'gain');
    signal.Powerz = W2dBm(signal.Powerz);
        signal.Powerz = signal.Powerz - signal.Powerz(1);
    pumps.Powerz  = W2dBm(pumps.Powerz);
        pumps.Powerz = pumps.Powerz - pumps.Powerz(1,:).*pumps.FWDchannels;
        pumps.Powerz = pumps.Powerz - pumps.Powerz(end,:).*pumps.BKWchannels;
    yMax = max([max(signal.Powerz),max(pumps.Powerz)]);
    if ~isfield(options,'ymin'), options.ymin = yMax - 10;end
    if ~isfield(options,'ymax'), options.ymax = yMax;end
    ystr = 'Gain (dB)';
    yL = [options.ymin, options.ymax];
elseif strcmpi(options.units,'dBm')
    signal.Powerz = W2dBm(signal.Powerz);
    pumps.Powerz  = W2dBm(pumps.Powerz);
    yMax = max([max(signal.Powerz),max(pumps.Powerz)]);
    if ~isfield(options,'ymin'), options.ymin = yMax - 10;end
    if ~isfield(options,'ymax'), options.ymax = yMax;end
    ystr = 'Gain (dB)';
    yL = [options.ymin, options.ymax];
    ystr = 'Power (dBm)';
else strcmpi(options.units,'linear')
    ystr = 'Power (W)';
    if ~isfield(options,'ymin'), options.ymin = 0;end
    if ~isfield(options,'ymax'), options.ymax = max([max(signal.Powerz),max(pumps.Powerz)]).*1.2;end
    yL = [options.ymin, options.ymax];
end
if max(yL) == 0 || max(yL) == -Inf, yL= [0,1];end

figure(sys.fignum);clf;
plot(sys.z,signal.Powerz,'-b','linewidth',3)
hold on; plot(sys.z,pumps.Powerz(:,pumps.FWDchannels),'-r','linewidth',3)
hold on; plot(sys.z,pumps.Powerz(:,pumps.BKWchannels),'--r','linewidth',3)
ylim(yL)
xlim([0,fiber.Lz]);grid on;
ylabel(ystr);
yyaxis('right');
plot(sys.z(2:end),fiber.N2z.*100,'--k')
ylabel('N_2 (%)');
set(gca,'Ycolor',[0,0,0])
ylim([0,100])
xlabel('Propagation (m)');
legend('Signal','P_{fw}','P_{bw}','n_2')
sys.fignum = sys.fignum+1;


end