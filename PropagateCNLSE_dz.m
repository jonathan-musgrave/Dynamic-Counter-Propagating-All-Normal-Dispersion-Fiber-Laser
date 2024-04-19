
   % update signal.h, signal.v, signal.H, signal.V, signal.PSD,
    % signal.PulseEng, signal.Peakpower for 1 step dz in m
function signal = PropagateCNLSE_dz(sys, signal, fiber, Gsat,  dz)
    % reassign variables for easy of reading
    dt = sys.dt;
    nt = sys.nt;
    frep = sys.frep;
    
    % BC((abs(signal.t)./1e-12)>20) = 0;
    %% Linear step
    h_temp = signal.h;
    v_temp = signal.v;
    
    %% Aply dispersion and gain to each polarization seperatley
    signal.H = IFFT(h_temp,nt,dt,frep).*exp(fiber.dispersion*(1j*dz)).*exp(Gsat/2.*dz);
    signal.V = IFFT(v_temp,nt,dt,frep).*exp(fiber.dispersion*(1j*dz)).*exp(Gsat/2.*dz);
    % Back to Time domain;
    h = FFT(signal.H,nt,dt,frep);
    v = FFT(signal.V,nt,dt,frep);
    
%     
%     h = fft(signal.H);
%     v = fft(signal.V);
    %% Nonlinear step using RK4
    fiber.NonLin = -1j*fiber.gamma;
    % h = fft(signal.H); v = ifft(signal.V);
    % RK4 Step for dz
    AA11=(abs(h).^2+2/3*abs(v).^2).*h+ 1/3*v.*v.*conj(h);
    kh1 = -fiber.NonLin*(AA11);% +(-1i)/(281.7*2*pi)*deriv1(AA11,dt));% ifft(1i*signal.w.*fft(AA11)));
    
    AA12=(abs(v).^2+2/3*abs(h).^2).*v + 1/3*h.*h.*conj(v);
    kv1 = -fiber.NonLin*(AA12);%+(-1i)/(281.7*2*pi)*deriv1(AA12,dt));% *ifft(1i*signal.w.*fft(AA12)));
   
    uh_half2 = h + kh1*dz/2;
    uv_half2 = v + kv1*dz/2;
    
    
    AA21=(abs(uh_half2).^2+2/3*abs(uv_half2).^2).*uh_half2 + 1/3*uv_half2.*uv_half2.*conj(uh_half2);
    kh2 = -fiber.NonLin*(AA21);% +(-1i)/(281.7*2*pi)*deriv1(AA21,dt));%*ifft(1i*signal.w.*fft(AA21)));
    
    AA22=(abs(uv_half2).^2+2/3*abs(uh_half2).^2).*uv_half2 + 1/3*uh_half2.*uh_half2.*conj(uv_half2);
    kv2 = -fiber.NonLin*(AA22);%+(-1i)/(281.7*2*pi)*deriv1(AA22,dt));%*ifft(1i*signal.w.*fft(AA22)));
    
    uh_half3 = h + kh2*dz/2;
    uv_half3 = v + kv2*dz/2;
    
    AA31=(abs(uh_half3).^2+2/3*abs(uv_half3).^2).*uh_half3 + 1/3*uv_half3.*uv_half3.*conj(uh_half3);
    kh3 = -fiber.NonLin*(AA31);%+(-1i)/(281.7*2*pi)*ifft(1i*signal.w.*fft(AA31)));
    
    AA32=(abs(uv_half3).^2+2/3*abs(uh_half3).^2).*uv_half3 + 1/3*uh_half3.*uh_half3.*conj(uv_half3);
    kv3 = -fiber.NonLin*(AA32);%+(-1i)/(281.7*2*pi)*deriv1(AA32,dt));%*ifft(1i*signal.w.*fft(AA32)));
   
    uh_full = h + kh3*dz;
    uv_full = v + kv3*dz;
    
    AA41=(abs(uh_full).^2+2/3*abs(uv_full).^2).*uh_full + 1/3*uv_full.*uv_full.*conj(uh_full);
    kh4 = -fiber.NonLin*(AA41);% +(-1i)/(281.7*2*pi)*deriv1(AA41,dt));%*ifft(1i*signal.w.*fft(AA41)));
    
    AA42=(abs(uv_full).^2+2/3*abs(uh_full).^2).*uv_full + 1/3*uh_full.*uh_full.*conj(uv_full);
    kv4 = -fiber.NonLin*(AA42);% +(-1i)/(281.7*2*pi)*deriv1(AA42,dt));%*ifft(1i*signal.w.*fft(AA42)));
    
    % Update signal.h and signal.v
    signal.h = (h + dz/6*(kh1+2*kh2+2*kh3+kh4));
    signal.v = (v + dz/6*(kv1+2*kv2+2*kv3+kv4));
    signal.H = IFFT(signal.h,nt,dt,frep);
    signal.V = IFFT(signal.v,nt,dt,frep);
    
    % Resave values
    signal.PSD = abs(signal.H).^2+abs(signal.V).^2; signal.Power = sum(signal.PSD);
    signal.PulseEng = sum(abs(signal.h).^2+abs(signal.v).^2).*dt; % Calculate the updated pulse energy for step dz;
    signal.PeakPower = max(abs(signal.h).^2+abs(signal.v).^2);
    
    
end