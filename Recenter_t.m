function [signal, dtShift, frep]  = Recenter_t(signal, sys, RTShiftCounter);

DT = sum(RTShiftCounter); % sums all the delta ts of each round trip
DT = 0;
uuOut = signal.uuOut; % Change this to whatever you want to recenter
uu    = signal.uu; % [2xnt];
t = sys.t;
com_t = calc_COM(t,abs(uuOut).^2);

[~,inx] = min(abs(com_t-t));

uu_shift = fftshift(circshift(uu,[0,-inx]),2); % Circshift COM to inx = 0  than fftshift to center at t = 0

frep = 1./(1./(sys.frep) - DT - com_t); % note minus sign because negative in time window is forward side of the pulse
dtShift = com_t;


% Resave all the values with the new shifted version;
signal.uu = uu_shift;
signal.uuOUt = fftshift(circshift(uuOut,-inx,2));


signal.h = signal.uu(1,:); signal.v = signal.uu(2,:); % Update CW signal; (H and V are updated automattically in Fiber_Prop_Lz)
    


end
